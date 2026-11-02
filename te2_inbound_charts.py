# -*- coding: utf-8 -*-
"""
T/E2 ADJUST — 一站式入站流水线
================================
整个作业由单个脚本完成：

  1. 抓取 README.html，解析最新的“逐日明细”页面，即年份最大的
     Inbound-<year>.html 链接（当统计窗口跨年时，回退到 Inbound-<year-1>.html）；
  2. 解析每日记录，取所需统计周期内的阿那曲唑 / 来曲唑剂量；
  3. 运行一室一级吸收 PK 模型（对每一剂做线性叠加），绘制最终的双轴图表。

输出：
  默认：<output.dir>/Inbound-<结束日期yyyymmdd>-<统计天数>.<默认扩展名>
  以及同目录下与图片同名的 .json 元数据文件。

用法
----
  python te2_inbound_charts.py [--config PATH]

所有可调参数均来自 TOML 配置文件：
  [run].end / [run].period / [run].window
  [output].dir / [output].path
  [site].proxy
默认配置文件名为“脚本名去掉扩展名 + .toml”，
例如脚本为 te2_inbound_charts.py 时，配置文件名是 te2_inbound_charts.toml。
命令行只保留 --config，用于在需要时显式指定配置文件路径。

说明
----
* 需要 Python 3.12 或更高版本。
* 网络：通过 urllib 抓取；是否使用代理由 [site].proxy 决定
  （为空时显式禁用代理，非空时作为 HTTP/HTTPS 代理使用）。
* 页面缓存在系统临时目录：<temp>/<site.cache_dir_name>/。
* 每日 08:00 的采样为给药前谷值，与既定规范一致。
"""
# region # <import>  引用扩展库
import sys
if sys.version_info < (3, 12):
    raise SystemExit("Python 3.12 或更高版本是必需的。")
import os, tomllib, re, tempfile, ssl, urllib.request as urlreq
import time, math, argparse, json, html as html_mod, numpy as np
from numpy.typing import NDArray
from datetime import datetime, timedelta, date
from scipy.interpolate import UnivariateSpline, make_interp_spline
import matplotlib, matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
from matplotlib.artist import Artist
# endregion
# region # <alias>   类型别名
type ListDT = list[datetime]
type LStr   = list[str]
type LInt   = list[int]
type LFloat = list[float]
type SStr   = set[str]
type NPFloats = NDArray[np.float64]
type DFloat = dict[str, float]
type DInt   = dict[str, int]
type DStr   = dict[str, str]
type DictParam = dict[str, str|int]
type DataRow   = dict[str, str|float]
type DataRows  = list[DataRow]
type DataPages = dict[int, DataRows]
type DataDose  = tuple[datetime, float|str]
type DataDoses = list[DataDose]
type DBDoses   = dict[str, DataDoses]
type DataAxes  = dict[str, float|NPFloats]
type DataImg   = dict[str, LStr|str]
type UrlProxy   = urlreq.OpenerDirector
type UrlReq     = urlreq.Request
type UrlHandler = urlreq.ProxyHandler
type ArgName    = argparse.Namespace
type ArgParser  = argparse.ArgumentParser
type TConfig = dict[str, dict]
type TMeta   = dict[str, bool|str|dict]
# endregion
# region # <runtime> 运行时配置（由 load_config 从 TOML 填充）
CFG: TConfig
_SCRIPT_BASENAME: str = os.path.basename(__file__)
_SCRIPT_STEM: str = os.path.splitext(_SCRIPT_BASENAME)[0]
DEFAULT_CONFIG_NAME: str = _SCRIPT_STEM + ".toml"

# V_per_F == 0 时必须同时提供的标定参数
PK_CALIB_KEYS: tuple[str, str, str] = ("calib_trough", "calib_dose", "calib_tau")

def log(*content: str) -> None:
    """统一日志输出（立即刷新）。"""
    print(*content, flush=True)
# endregion
# region # 0. 配置加载
def default_config_path() -> str:
    """默认配置文件：与本脚本同目录的 <脚本名去掉扩展名>.toml。"""
    here: str = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, DEFAULT_CONFIG_NAME)

def _validate_pk_section(pk: DFloat, drug: str, path: str) -> None:
    """校验单个 [pk.<drug>] 段。
    V_per_F > 0 时直接使用该值；V_per_F == 0（或缺失）时表示需要按谷值
    反算，此时 calib_trough / calib_dose / calib_tau 三项都必须提供。
    """
    if not isinstance(pk, dict):
        raise SystemExit("配置文件 [pk.%s] 不是有效表：%s" % (drug, path))

    for k in ("t_half_h", "ka"):
        if k not in pk:
            raise SystemExit("配置文件 [pk.%s] 缺少必需参数 %s：%s" % (drug, k, path))

    v: float = pk.get("V_per_F", 0.0) or 0.0
    if v > 0.0: return # 直接给定 V/F，不需要 calib_*

    missing: LStr = [k for k in PK_CALIB_KEYS if k not in pk]
    if missing: raise SystemExit()
    if missing: raise SystemExit(
        "配置文件 [pk.%s] 的 V_per_F 为 0（表示按谷值反算），"
        "但缺少标定参数：%s：%s"
        % (drug, ", ".join(missing), path)
    )

def load_config(path: str) -> TConfig:
    """从 TOML 文件读取全部配置，并做基本校验。"""
    if not os.path.isfile(path):
        raise SystemExit("配置文件不存在：%s" % path)
    with open(path, "rb") as f:
        try:
            data: TConfig = tomllib.load(f)
        except tomllib.TOMLDecodeError as e:
            raise SystemExit("配置文件解析失败：%s (%s)" % (path, e))

    required: tuple[str, ...] = ("app", "site", "run", "output", "render",
                                 "colors", "fonts", "image", "period", "pk")
    for s in required:
        if s not in data:
            raise SystemExit("配置文件缺少段 [%s]：%s" % (s, path))

    for drug in ("anastrozole", "letrozole"):
        if drug not in data["pk"]:
            raise SystemExit("配置文件缺少 [pk.%s] 段：%s" % (drug, path))
        _validate_pk_section(data["pk"][drug], drug, path)
    return data

def parse_period(value: str|int) -> int:
    """解析 [run].period（大小写不敏感的别名，或整数天数）。"""
    s: str = str(value).strip().lower()
    p: DInt = CFG["period"]
    if s in ("monthly", "quarter", "half", "year"): return int(p[s])
    if re.fullmatch(r"\d+", s): return int(s)
    raise SystemExit("配置错误：[run].period 必须为 monthly、quarter、half、"
                     "year 或整数天数，实际为 %r" % value)

def load_run_params() -> DictParam:
    """从 [run] / [output] / [site] 读取并校验全部运行参数。"""
    run:  DictParam = CFG["run"]
    out:  DStr   = CFG["output"]
    site: DStr   = CFG["site"]
    p:    DInt   = CFG["period"]

    period: int = parse_period(run["period"])
    window: int = int(run["window"])

    if window < p["min_window_days"]:
        raise SystemExit("配置错误：[run].window 必须 >= %d"
                         % int(p["min_window_days"]))
    if period < window:
        raise SystemExit("配置错误：[run].period 天数必须 >= [run].window")
    if period > int(p["max_days"]):
        raise SystemExit("配置错误：[run].period 天数必须 <= %d"
                         % int(p["max_days"]))

    return {
        "end": str(run["end"]),
        "period": period,
        "window": window,
        "output_dir": out["dir"],
        "output_path": out["path"],
        "proxy": str(site.get("proxy", "") or ""),
    }
# endregion
# region # 1. 抓取
def _cache_dir() -> str:
    """缓存目录：系统临时目录下的 cache_dir_name。"""
    return os.path.join(tempfile.gettempdir(), CFG["site"]["cache_dir_name"])

def _normalize_proxy(proxy: str) -> str:
    """规范化代理地址；空字符串表示不使用代理。"""
    p: str = (proxy or "").strip()
    if not p: return ""
    if "://" not in p: p = "http://" + p
    return p

def _opener(proxy: str = "") -> UrlProxy:
    """构造 urllib opener；proxy 为空时显式禁用代理。"""
    ctx: ssl.SSLContext = ssl.create_default_context()
    p: str = _normalize_proxy(proxy)
    proxy_handler: UrlHandler
    if p: proxy_handler = urlreq.ProxyHandler({"http": p, "https": p})
    else: proxy_handler = urlreq.ProxyHandler({}) # 显式禁用环境变量中的代理
    return urlreq.build_opener(proxy_handler, urlreq.HTTPSHandler(context=ctx))

def fetch(name: str, timeout: int = 60, proxy: str = "") -> tuple[str|None, str]:
    """从站点抓取 <name>，缓存到系统临时目录；网络失败时回退到缓存。"""
    cache: str = _cache_dir()
    path: str = os.path.join(cache, name)
    url: str = CFG["site"]["base"] + name
    t0: float = time.time()
    try:
        req: UrlReq = urlreq.Request(url, headers={"User-Agent": CFG["site"]["user_agent"]})
        with _opener(proxy).open(req, timeout=timeout) as r: raw = r.read()
        os.makedirs(cache, exist_ok=True)
        with open(path, "wb") as f: f.write(raw)
        log("      GET %-24s %6d bytes  %.2fs" % (name, len(raw), time.time() - t0))
        return raw.decode("utf-8", "replace"), "network"
    except Exception as e:
        log("      GET %-24s FAILED (%s: %s)" % (name, type(e).__name__, e))
    if os.path.exists(path):
        txt: str = open(path, encoding="utf-8", errors="replace").read()
        log("      CACHE %-22s %6d bytes" % (name, os.path.getsize(path)))
        return txt, "cache"
    return None, "missing"

def newest_inbound_year(readme_html) -> tuple[int, LInt]:
    """从 README 解析年份最大的 Inbound-<year>.html 链接。"""
    years: LInt = sorted({int(y) for y in re.findall(r'Inbound-(\d{4})\.html', readme_html)})
    if not years: raise SystemExit("README 中未找到任何 Inbound-<year>.html 链接")
    return years[-1], years
# endregion
# region # 2. 解析每日记录
def cell_text(td: str) -> str:
    """提取表格单元格的纯文本。"""
    td = re.sub(r"(?is)<br\s*/?>", " / ", td)
    td = re.sub(r"(?s)<[^>]+>", " ", td)
    td = html_mod.unescape(td).replace("\u00a0", " ")
    return re.sub(r"\s+", " ", td).strip()

def num(s: str) -> float | None:
    """把单元格文本转换为数字，缺失值返回 None。"""
    s = (s or "").strip()
    if s in ("", "-", "--"): return None
    return float(s) if re.match(r"^-?\d+(?:\.\d+)?$", s) else None

def parse_inbound(page_html: str) -> DataRows:
    """-> 按日期升序排列的 [{date, anastrozole_mg, letrozole_mg, ...}]"""
    rows: DataRows = []
    seen: SStr = set()
    for tr in re.findall(r"(?is)<tr>(.*?)</tr>", page_html):
        tds: LStr = re.findall(r"(?is)<td[^>]*>(.*?)</td>", tr)
        if len(tds) < 6:
            continue
        v: LStr = [cell_text(c) for c in tds]
        if not re.match(r"^\d{4}-\d{2}-\d{2}$", v[0]) or v[0] in seen:
            continue
        seen.add(v[0])
        rows.append({"date": v[0] or "",
                     "anastrozole_mg": num(v[2]) or 0,
                     "letrozole_mg": num(v[3]) or 0,
                     "testosterone_mg": (num(v[4]) if len(v) > 4 else 4) or 0,
                     "hcg_IU": (num(v[5]) if len(v) > 5 else 5) or 0})
    rows.sort(key=lambda r: r["date"])
    return rows
# endregion
# region # 3. PK 模型（一室、一级吸收、线性叠加）
class Drug:
    """一室一级吸收模型。"""
    def __init__(self, key: str) -> None:
        p: DFloat = dict(CFG["pk"][key])
        self.name: str = key
        self.ke: float = math.log(2.0) / p["t_half_h"]
        self.ka: float = p["ka"]
        self.V: float|None = float(p.get("V_per_F", 0.0) or 0.0) or None
        if self.V is None:
            # V_per_F == 0：按说明书谷值反算 V/F
            # （标定参数的完整性在 load_config 中已校验）
            R: float = math.exp(-self.ke * p["calib_tau"])
            shape: float = (self.ka / (self.ka - self.ke)) * \
                    (math.exp(-self.ke * p["calib_tau"]) -
                     math.exp(-self.ka * p["calib_tau"]))
            self.V = 1000.0 * p["calib_dose"] * shape / (1.0 - R) / p["calib_trough"]

    def conc(self, dt_h: float, dose_mg: float) -> float:
        """单剂在 dt_h 小时后的血药浓度贡献。"""
        if dose_mg is None or dose_mg <= 0 or dt_h < 0: return 0.0
        c: float = (dose_mg * self.ka) / (self.V * (self.ka - self.ke))
        return 1000.0 * c * (math.exp(-self.ke * dt_h) - math.exp(-self.ka * dt_h))

    def concentration(self, when: datetime, doses: DataDoses) -> float:
        """对回溯窗口内的每一剂做线性叠加。
        剂量回溯窗口很关键：超过约 20 个半衰期的一剂仍会贡献 ~1e-5 ng/mL，
        足以在第 4 位小数上翻转结果（进而影响抗锯齿像素），
        与截断式仿真产生可见差异。
        """
        total: float = 0.0
        lookback: int = CFG["render"]["dose_lookback_h"]
        for t, mg in doses:
            dt_h: float = (when - t).total_seconds() / 3600.0
            if dt_h < 0 or dt_h > lookback:
                continue
            total += self.conc(dt_h, float(mg))
        return total

def build_doses(rows: DataRows, field: str) -> DataDoses:
    """从每日记录中提取给药时间点与剂量。"""
    hour: int = CFG["render"]["sample_hour"]
    return [(datetime.strptime(str(r["date"]), "%Y-%m-%d") + timedelta(hours=hour), r[field])
            for r in rows if r.get(field)]

def daily_series(drug: Drug, doses: DataDoses,
                 start: datetime, end: datetime) -> NPFloats:
    """按每天一次、在 sample_hour 采样，首尾都包含。"""
    step: float = CFG["render"]["sample_step_h"]
    out: LFloat = []
    t: datetime = start
    while t <= end:
        out.append(drug.concentration(t, doses))
        t += timedelta(hours=step)
    return np.array(out)

def canon(series: NPFloats, decimals: int|None = None) -> NPFloats:
    """按 CSV 存储的方式保留小数：格式化成字符串再解析回来，
    保证结果与 CSV 往返后完全一致（np.round 可能偏移 1 ULP 并影响抗锯齿像素）。"""
    if decimals is None: decimals = CFG["render"]["conc_decimals"]
    fmt: str = "%." + str(decimals) + "f"
    return np.array([float(fmt % v) for v in series])

def rolling_vol(x: NPFloats, window: int) -> NPFloats:
    """滚动窗口标准差。"""
    out: NPFloats = np.full(len(x), np.nan)
    for i in range(window - 1, len(x)):
        out[i] = x[i - window + 1: i + 1].std(ddof=1)
    return out

def fit_bar_curve(x: NPFloats, y: NPFloats) -> tuple[UnivariateSpline, DataRow]:
    """按差分的经验标准差选取平滑参数，拟合色带顶部曲线。"""
    sigma: float = np.std(np.diff(y), ddof=1) / np.sqrt(2.0)
    s: float = len(x) * sigma ** 2
    spl: UnivariateSpline = UnivariateSpline(x, y, k=3, s=s)
    return spl, {"sigma_est": round(sigma, 6), "s": round(s, 6),
                 "residual_rms": round(np.sqrt(
                     np.mean((y - np.clip(spl(x), 0, None)) ** 2)), 6),
                 "raw_mean": round(float(y.mean()), 6),
                 "fitted_mean": round(float(np.clip(spl(x), 0, None).mean()), 6)}
# endregion
# region # 4. 渲染
def setup_fonts() -> str:
    """按配置设置等宽字体，返回实际使用的字体名。"""
    from matplotlib.font_manager import fontManager
    prefer: LStr = CFG["fonts"]["monospace"]
    have: SStr = {f.name for f in fontManager.ttflist}
    picked: LStr = [f for f in prefer if f in have]
    matplotlib.rcParams["font.family"] = "monospace"
    matplotlib.rcParams["font.monospace"] = picked + ["DejaVu Sans Mono"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    return picked[0] if picked else "DejaVu Sans Mono"

def render(wd: ListDT, a_win: NPFloats, l_win: NPFloats, b_ana: NPFloats,
           b_let: NPFloats, xf: NPFloats, a_smooth: NPFloats, l_smooth: NPFloats,
           left_ticks: float|NPFloats, right_ticks: float|NPFloats,
           ymax: float, right_top: float, xlabel: str) -> Figure:
    """生成双轴图表。"""
    c: DataRow = CFG["colors"]
    app: DStr = CFG["app"]
    n_int: int = CFG["render"]["n_intervals"]

    fig: Figure; ax: Axes; ax2: Axes
    fig, ax = plt.subplots(figsize=(16.0, 9.0), dpi=110)
    fig.patch.set_facecolor("white")
    ax2 = ax.twinx()
    ax2.set_zorder(ax.get_zorder() - 1)
    ax.patch.set_visible(False)

    xs: NPFloats = np.linspace(xf[0], xf[-1], len(b_ana))
    ax2.fill_between(xs, 0.0, b_ana, color=c["ana_bar"], linewidth=0, zorder=2)
    ax2.fill_between(xs, b_ana, b_ana + b_let, color=c["let_bar"], linewidth=0, zorder=2)
    ax2.set_ylim(0, right_top)

    ax.plot(xf, a_smooth, color=c["ana_line"], linewidth=2.6,
            solid_capstyle="round", label=app["label_anastrozole"], zorder=5)
    ax.plot(xf, l_smooth, color=c["let_line"], linewidth=2.6,
            solid_capstyle="round", label=app["label_letrozole"], zorder=5)
    x: NPFloats = np.arange(len(wd), dtype=float)
    ax.scatter(x, a_win, s=7, color=c["ana_line"], alpha=0.28, zorder=4)
    ax.scatter(x, l_win, s=7, color=c["let_line"], alpha=0.28, zorder=4)

    ax.set_ylim(0, ymax)
    ax.set_xlim(-0.6, len(wd) - 0.4)
    ax.set_yticks(left_ticks)
    ax2.set_yticks(right_ticks)
    # 右轴标签小数位随步长自适应，避免小范围图表退化为「0 / 1」
    r_step: float = right_top / n_int
    r_dec: int = 0 if r_step >= 1.0 else (2 if r_step < 0.1 else 1)
    ticks_arr: NPFloats = np.atleast_1d(right_ticks)
    ax2.set_yticklabels([f"{v:.{r_dec}f}" for v in ticks_arr])

    ax.set_ylabel(app["ylabel_left"], fontsize=12.5, color=c["axis_text_left"])
    ax2.set_ylabel(app["ylabel_right"], fontsize=12.5, color=c["axis_text_right"])
    ax.tick_params(axis="y", colors=c["axis_text_left"], labelsize=10)
    ax2.tick_params(axis="y", colors=c["axis_text_right"], labelsize=10)

    tick_idx: LInt = [i for i, d in enumerate[datetime](wd) if d.day in (1, 16)]
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([wd[i].strftime("%m-%d") for i in tick_idx], fontsize=10)
    ax.set_xlabel(xlabel, fontsize=11.5)
    for i in tick_idx: ax.axvline(i, color=c["grid"], linewidth=c["grid_lw"], zorder=3)
    for v in left_ticks[1:]: ax.axhline(v, color=c["grid"], linewidth=c["grid_lw"], zorder=3)
    ax.set_axisbelow(False)

    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax.spines["left"].set_color(str(c["spine"]))
    ax2.spines["right"].set_color(str(c["spine"]))

    ax.set_title(app["title"], fontsize=16, fontweight="bold", pad=14)
    h1: list[Artist]; l1: LStr
    h1, l1 = ax.get_legend_handles_labels()
    ax.legend(h1, l1, loc="upper left", fontsize=10.5, ncol=2, framealpha=0.0)

    fig.tight_layout()
    return fig
# endregion
# region # 5. 输出路径
def resolve_output_path(output_dir: str, output_path: str, default_name: str) -> str:
    """解析 [output].dir 与 [output].path，兼容 Windows 与 Linux。
    - output_dir 可为相对或绝对路径，统一规范化为绝对路径；
    - output_path 为空 -> output_dir/default_name；
    - output_path 不含目录部分 -> 放入 output_dir；
    - output_path 相对且含目录 -> 相对于当前工作目录；
    - output_path 为绝对路径 -> 原样使用；
    - 无扩展名 -> 自动补默认扩展名；
    - 扩展名不在支持列表 -> 抛出 ValueError。
    """
    output_dir = os.path.abspath(os.path.expanduser(output_dir or "."))
    if output_path is None or str(output_path).strip() == "":
        return os.path.join(output_dir, default_name)

    p: str = os.path.expanduser(str(output_path).strip())
    resolved: str
    if os.path.isabs(p): resolved = os.path.normpath(p)
    elif os.path.dirname(p): resolved = os.path.abspath(p)
    else: resolved = os.path.join(output_dir, p)

    img: DataImg = CFG["image"]
    supported: tuple[str, ...] = tuple[str, ...](e.lower() for e in img["supported_exts"])
    default_ext: str = img["default_ext"].lower()

    ext: str = os.path.splitext(resolved)[1].lower()
    if not ext: resolved += default_ext
    elif ext not in supported:
        raise ValueError(
            "[output].path 的扩展名 %r 不受支持；可用：%s"
            % (ext, ", ".join(supported))
        )
    return resolved
# endregion
# region # 6. 流水线各阶段
def resolve_page_year(proxy: str) -> tuple[int, str]:
    """抓取 README.html，返回（最新年份, README 来源）。"""
    readme: str|None; src_r: str
    readme, src_r = fetch("README.html", proxy=proxy)
    if readme is None:
        raise SystemExit("README.html 不可用（网络与缓存均失败）")
    year: int; all_years: LInt
    year, all_years = newest_inbound_year(readme)
    log("      README source=%s  Inbound years seen=%s  -> newest = %d"
        % (src_r, all_years, year))
    return year, src_r

def fetch_main_page(year: int, proxy: str) -> DataPages:
    """抓取最新的 Inbound-<year>.html，返回 {year: rows}。"""
    target: str = "Inbound-%d.html" % year
    txt: str|None; src: str
    txt, src = fetch(target, proxy=proxy)
    if txt is None: raise SystemExit("%s 不可用" % target)
    pages: DataPages = {year: parse_inbound(txt)}
    log("      %s (%s): %d daily rows  %s .. %s"
        % (target, src, len(pages[year]),
           pages[year][0]["date"], pages[year][-1]["date"]))
    return pages

def find_last_recorded(rows: DataRows) -> str:
    """最后一条包含剂量的记录日期。"""
    return str(next(r["date"] for r in reversed(rows)
                if r["anastrozole_mg"] or r["letrozole_mg"]))

def resolve_end_date(last_recorded: str, end_arg: str) -> date:
    """根据 [run].end 决定窗口结束日。"""
    if end_arg == "auto": return datetime.strptime(last_recorded, "%Y-%m-%d").date()
    return datetime.strptime(end_arg, "%Y-%m-%d").date()

def compute_window(end_date: date, period_days: int, window_days: int) \
    -> tuple[datetime, datetime, datetime]:
    """根据统计周期计算 (win_start, win_end, sim_start)。"""
    hour: int = CFG["render"]["sample_hour"]
    win_end: datetime = datetime.combine(end_date, datetime.min.time()) + timedelta(hours=hour)
    win_start: datetime = win_end - timedelta(days=period_days - 1)
    # 额外天数用于覆盖滚动窗口 + PK 剂量回溯
    sim_start: datetime = win_start - timedelta(days=window_days - 1 + 40)
    return win_start, win_end, sim_start

def maybe_extend_previous_year(rows: DataRows, year: int, sim_start: datetime,
                               proxy: str, pages: DataPages) -> DataRows:
    """仿真起点跨年时，将上一年的记录前置到 rows 之前；pages 就地更新。"""
    if sim_start.year >= year: return rows
    prev: str = "Inbound-%d.html" % (year - 1)
    t2: str|None; s2: str
    t2, s2 = fetch(prev, proxy=proxy)
    if t2 is None:
        log("      %s 不可用——%d-01-01 之前的剂量未知，"
            "最初数周可能偏低" % (prev, year))
        return rows
    pages[year - 1] = parse_inbound(t2)
    merged: DataRows = pages[year - 1] + rows
    log("      %s (%s): %d daily rows %s .. %s"
        % (prev, s2, len(pages[year - 1]),
           pages[year - 1][0]["date"], pages[year - 1][-1]["date"]))
    merged.sort(key=lambda r: r["date"])
    return merged

def build_dose_table(rows: DataRows) -> DBDoses:
    """构建两种药物的给药列表并打印概览。"""
    doses: DBDoses = {k: build_doses(rows, k) for k in ("anastrozole_mg", "letrozole_mg")}
    for k, v in doses.items():
        log("      %-16s doses=%3d  total=%7.1f mg  first=%s"
            % (k, len(v), sum(float(m) for _, m in v), v[0][0].date() if v else "-"))
    return doses

def simulate_pk(doses: DBDoses, sim_start: datetime, win_end: datetime) \
        -> tuple[Drug, Drug, NPFloats, NPFloats]:
    """运行 PK 仿真，返回两个药物对象以及完整浓度序列。"""
    ana: Drug = Drug("anastrozole")
    let: Drug = Drug("letrozole")
    log("      anastrozole t1/2=%.1fh, ka=%.2f, V/F=%.1f L" %
        (math.log(2) / ana.ke, ana.ka, ( 1 if ana.V is None else ana.V)))
    log("      letrozole   t1/2=%.1fh, ka=%.2f, V/F=%.1f L" %
        (math.log(2) / let.ke, let.ka, ( 1 if let.V is None else let.V)))
    a_full: NPFloats = canon(daily_series(ana, doses["anastrozole_mg"], sim_start, win_end))
    l_full: NPFloats = canon(daily_series(let, doses["letrozole_mg"], sim_start, win_end))
    log("      simulated points=%d" % len(a_full))
    return ana, let, a_full, l_full

def compute_bands(a_full: NPFloats, l_full: NPFloats, off: int, n: int, window: int) \
        -> tuple[NPFloats, NPFloats, NPFloats, NPFloats, DataRow, DataRow]:
    """计算滚动波动率、色带顶部平滑曲线和拟合诊断。"""
    v_ana: NPFloats = np.nan_to_num(rolling_vol(a_full, window))[off:]
    v_let: NPFloats = np.nan_to_num(rolling_vol(l_full, window))[off:]
    x: NPFloats = np.arange(n, dtype=float)
    spl_a: UnivariateSpline; info_a: DataRow
    spl_a, info_a = fit_bar_curve(x, v_ana)
    spl_l: UnivariateSpline; info_l: DataRow
    spl_l, info_l = fit_bar_curve(x, v_let)
    sub: int = CFG["render"]["sub_per_day"]
    xs_sub: NPFloats = np.linspace(x[0], x[-1], n * sub)
    top_a: NPFloats = np.clip(spl_a(xs_sub), 0.0, None)
    top_l: NPFloats = np.clip(spl_l(xs_sub), 0.0, None)
    return v_ana, v_let, top_a, top_l, info_a, info_l

def compute_axes(a_win: NPFloats, l_win: NPFloats, top_a: NPFloats, top_l: NPFloats) -> DataAxes:
    """计算两个 y 轴共用的刻度与上限几何。"""
    n_int: int = CFG["render"]["n_intervals"]
    mean_h: float = (top_a + top_l).mean()
    right_top: float = mean_h / 0.25
    data_max: float = max(a_win.max(), l_win.max())
    step_l: float = math.ceil(data_max / n_int)
    ymax: float = step_l * n_int
    left_ticks: NPFloats = np.linspace(0.0, ymax, n_int + 1)
    right_ticks: NPFloats = np.linspace(0.0, right_top, n_int + 1)
    return {"mean_h": mean_h, "right_top": right_top, "step_l": step_l,
            "ymax": ymax, "left_ticks": left_ticks, "right_ticks": right_ticks}

def render_chart(wd: ListDT, a_win: NPFloats, l_win: NPFloats,
                 top_a: NPFloats, top_l: NPFloats, axes: DataAxes) -> Figure:
    """对真实采样做三次样条加密后调用 render。"""
    n: int = len(wd)
    sub: int = CFG["render"]["sub_per_day"]
    x: NPFloats = np.arange(n, dtype=float)
    xf: NPFloats = np.linspace(x[0], x[-1], n * sub)
    a_smooth: NPFloats = np.clip(make_interp_spline(x, a_win, k=3)(xf), 0.0, None)
    l_smooth: NPFloats = np.clip(make_interp_spline(x, l_win, k=3)(xf), 0.0, None)
    return render(wd, a_win, l_win, top_a, top_l, xf, a_smooth, l_smooth,
                  axes["left_ticks"], axes["right_ticks"], float(axes["ymax"]),
                  float(axes["right_top"]), CFG["app"]["x_label"])

def save_figure(fig: Figure, out_path: str) -> None:
    """保存图表；输出格式由扩展名决定。"""
    fig.savefig(out_path, facecolor="white")
    plt.close(fig)
    log("      wrote %s (%d bytes)" % (out_path, os.path.getsize(out_path)))

def metadata_path(out_path: str) -> str:
    """与图片同名、同目录，扩展名替换为 .json。"""
    return os.path.splitext(out_path)[0] + ".json"

def write_metadata(out_path: str, wd: ListDT, last_recorded: str, src_r: str, pages: DataPages,
                   window: int, a_win: NPFloats, l_win: NPFloats, axes: DataAxes,
                   v_ana: NPFloats, v_let: NPFloats, top_a: NPFloats, top_l: NPFloats,
                   info_a: DataRow, info_l: DataRow, ana: Drug, let: Drug):
    """写出与图片同名的 JSON 元数据。"""
    app: DStr = CFG["app"]
    c: DataRow = CFG["colors"]
    n_int: int = CFG["render"]["n_intervals"]
    mean_h: float = float(axes["mean_h"])

    meta: TMeta = {
        "file": os.path.basename(out_path),
        "generated_by": _SCRIPT_BASENAME,
        "source": {"readme": src_r, "inbound_pages": sorted(pages.keys()),
                   "window_end_rule": "'auto' = 最后记录给药日 (%s)" % last_recorded},
        "window": {"start": wd[0].strftime("%Y-%m-%d"),
                   "end": wd[-1].strftime("%Y-%m-%d"),
                   "days": len(wd),
                   "sampling": "每天 08:00 采样一次（给药前谷值）"},
        "title": app["title"], "x_label": app["x_label"],
        "left_axis": {"label": app["ylabel_left"], "max": round(float(axes["ymax"]), 3),
                      "line": "穿过每个真实日采样点的三次插值样条",
                      "colors": {"anastrozole": c["ana_line"], "letrozole": c["let_line"]},
                      "peaks": {"anastrozole": round(float(a_win.max()), 3),
                                "letrozole": round(float(l_win.max()), 3)}},
        "right_axis": {"label": app["ylabel_right"], "max": round(float(axes["right_top"]), 4),
                       "mean_at_fraction": round(mean_h / axes["right_top"], 4),
                       "band_colors": {"anastrozole": c["ana_bar"], "letrozole": c["let_bar"]},
                       "bar_top": "用平滑拟合曲线替代旧的水平条顶段",
                       "moving_average_applied": False, "legend": False},
        "volatility": {"window_days": window, "statistic": "标准差",
                       "statistic_key": "std", "unit": "ng/mL",
                       "raw_mean_total": round(float((v_ana + v_let).mean()), 4),
                       "plotted_mean_total": round(mean_h, 4),
                       "raw_mean_abs_day_change": round(
                           float(np.abs(np.diff(v_ana + v_let)).mean()), 4),
                       "plotted_mean_abs_day_change": round(
                           float(np.abs(np.diff(top_a + top_l)).mean()), 4),
                       "fit_anastrozole": info_a, "fit_letrozole": info_l},
        "grid": {"vertical": "每月 1 日与 16 日",
                 "horizontal": "两个 y 轴、间隔数相同",
                 "color": c["grid"], "linewidth": c["grid_lw"], "zorder": 3},
        "axis_intervals": {"intervals_per_axis": n_int,
                           "left_step": round(float(axes["step_l"]), 4),
                           "left_max": round(float(axes["ymax"]), 4),
                           "right_step": round(float(axes["right_top"]) / n_int, 4),
                           "right_max": round(float(axes["right_top"]), 4)},
        "mean_reference_line": False,
        "stage_extrema": {"drawn": False, "note": "最终图表中不绘制极值标记"},
        "pk": {"anastrozole": {"t_half_h": CFG["pk"]["anastrozole"]["t_half_h"],
                               "ka": CFG["pk"]["anastrozole"]["ka"],
                               "V_per_F_L": round(ana.V or 0, 2),
                               "calibration": "稳态谷值 25.7 ng/mL @1 mg/日"},
               "letrozole": {"t_half_h": CFG["pk"]["letrozole"]["t_half_h"],
                             "ka": CFG["pk"]["letrozole"]["ka"],
                             "V_per_F_L": round(let.V or 0, 2)}},
    }
    meta_path: str = metadata_path(out_path)
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    log("      -> %s" % meta_path)
# endregion
# region # 7. 入口
def parse_args() -> ArgName:
    """命令行只保留 --config；其余参数全部来自 TOML。"""
    ap: ArgParser = argparse.ArgumentParser(
        description="TE2 入站图表流水线（所有参数见 TOML 配置文件）")
    ap.add_argument("--config", default=default_config_path(),
                    help="TOML 配置文件路径；默认与本脚本同目录的 %s"
                         % DEFAULT_CONFIG_NAME)
    return ap.parse_args()

def te2_inbound_charts_main():
    global CFG

    args: ArgName = parse_args()
    CFG = load_config(args.config)
    params: DictParam = load_run_params()
    period_days: int = int(params["period"])
    window: int = int(params["window"])
    proxy: str = str(params["proxy"])
    # 阶段 1：README + Inbound-<year>.html
    log("[1/6] resolving the newest daily-record page from README.html ...")
    year: int; src_r: str
    year, src_r = resolve_page_year(proxy)

    log("[2/6] fetching and parsing daily records ...")
    pages: DataPages = fetch_main_page(year, proxy)
    rows: DataRows = pages[year]

    # 阶段 2：窗口 + 跨年补齐
    last_recorded: str = find_last_recorded(rows)

    end_date: date = resolve_end_date(last_recorded, str(params["end"]))
    log("      last recorded dosing day = %s  -> window end = %s" % (last_recorded, end_date))

    log("[3/6] extending the window backwards ...")
    win_start: datetime; win_end: datetime; sim_start: datetime
    win_start, win_end, sim_start = compute_window(end_date, period_days, window)
    log("      window    %s .. %s  (%d days)"
        % (win_start.date(), win_end.date(), period_days))
    log("      simulated %s .. %s  (extra %d days for the rolling window and dose lookback)"
        % (sim_start.date(), win_end.date(), (win_start - sim_start).days))

    rows: DataRows = maybe_extend_previous_year(rows, year, sim_start, proxy, pages)
    doses: DBDoses = build_dose_table(rows)

    # 阶段 3：PK 仿真
    log("[4/6] PK simulation (1 point/day @%02d:00, pre-dose trough) ..."
        % CFG["render"]["sample_hour"])
    ana: Drug; let: Drug; a_full: NPFloats; l_full: NPFloats
    ana, let, a_full, l_full = simulate_pk(doses, sim_start, win_end)

    off: int = (win_start - sim_start).days
    a_win: NPFloats; l_win: NPFloats
    a_win, l_win = a_full[off:], l_full[off:]
    n: int = len(a_win)
    log("      window points=%d  peak ana=%.2f ng/mL  peak let=%.2f ng/mL"
        % (n, a_win.max(), l_win.max()))

    # 阶段 4：波动率 + 色带顶部
    log("[5/6] rolling %d-day volatility + smooth band top ..." % window)
    v_ana: NPFloats; v_let: NPFloats; top_a: NPFloats
    top_l: NPFloats; info_a: DataRow; info_l: DataRow
    v_ana, v_let, top_a, top_l, info_a, info_l = compute_bands(
        a_full, l_full, off, n, window)

    log("      fitted bar top |day-to-day|: raw=%.4f -> plotted=%.4f"
        % (np.abs(np.diff(v_ana + v_let)).mean(),
           np.abs(np.diff(top_a + top_l)).mean()))

    axes: DataAxes = compute_axes(a_win, l_win, top_a, top_l)
    log("      mean=%.4f -> right axis 0..%.4f (%.1f%%)  left axis 0..%.0f step=%.0f"
        % (float(axes["mean_h"]), float(axes["right_top"]),
           100 * axes["mean_h"] / axes["right_top"],
           float(axes["ymax"]), float(axes["step_l"])))

    # 阶段 5：渲染与输出
    log("[6/6] rendering ...")
    font: str = setup_fonts()
    log("      monospace font -> %s" % font)

    wd: ListDT = [win_start + timedelta(days=i) for i in range(n)]
    fig: Figure = render_chart(wd, a_win, l_win, top_a, top_l, axes)

    default_ext: str = CFG["image"]["default_ext"]
    default_name: str = "Inbound-%s-%d%s" % (wd[-1].strftime("%Y%m%d"), period_days, default_ext)
    try:
        out_path: str = resolve_output_path(str(params["output_dir"]),
                                            str(params["output_path"]), default_name)
    except ValueError as e:
        raise SystemExit("错误：%s" % e)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    save_figure(fig, out_path)

    write_metadata(out_path, wd, last_recorded, src_r, pages, window,
                   a_win, l_win, axes, v_ana, v_let, top_a, top_l,
                   info_a, info_l, ana, let)
    log("DONE")
# endregion
if __name__ == "__main__":
    te2_inbound_charts_main()

