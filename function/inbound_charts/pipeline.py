import os, math, json, numpy, matplotlib.pyplot
from datetime import date, datetime, timedelta
from scipy.interpolate import UnivariateSpline, make_interp_spline
from matplotlib.figure import Figure
from type_alias import LInt, LDT, DStr, TMeta, DataRow, DataRows, DataPages, DBDoses, NPFloats, DataAxes
from fetch import fetch, newest_inbound_year
from pk_model import Drug, build_doses, canon, daily_series, rolling_vol, fit_bar_curve
from parse_records import parse_inbound
from render import render
from runtime import log
import runtime

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
    hour: int = runtime.CFG["render"]["sample_hour"]
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
    v_ana: NPFloats = numpy.nan_to_num(rolling_vol(a_full, window))[off:]
    v_let: NPFloats = numpy.nan_to_num(rolling_vol(l_full, window))[off:]
    x: NPFloats = numpy.arange(n, dtype=float)
    spl_a: UnivariateSpline; info_a: DataRow
    spl_a, info_a = fit_bar_curve(x, v_ana)
    spl_l: UnivariateSpline; info_l: DataRow
    spl_l, info_l = fit_bar_curve(x, v_let)
    sub: int = runtime.CFG["render"]["sub_per_day"]
    xs_sub: NPFloats = numpy.linspace(x[0], x[-1], n * sub)
    top_a: NPFloats = numpy.clip(spl_a(xs_sub), 0.0, None)
    top_l: NPFloats = numpy.clip(spl_l(xs_sub), 0.0, None)
    return v_ana, v_let, top_a, top_l, info_a, info_l

def compute_axes(a_win: NPFloats, l_win: NPFloats, top_a: NPFloats, top_l: NPFloats) -> DataAxes:
    """计算两个 y 轴共用的刻度与上限几何。"""
    n_int: int = runtime.CFG["render"]["n_intervals"]
    mean_h: float = (top_a + top_l).mean()
    right_top: float = mean_h / 0.25
    data_max: float = max(a_win.max(), l_win.max())
    step_l: float = math.ceil(data_max / n_int)
    ymax: float = step_l * n_int
    left_ticks: NPFloats = numpy.linspace(0.0, ymax, n_int + 1)
    right_ticks: NPFloats = numpy.linspace(0.0, right_top, n_int + 1)
    return {"mean_h": mean_h, "right_top": right_top, "step_l": step_l,
            "ymax": ymax, "left_ticks": left_ticks, "right_ticks": right_ticks}

def render_chart(wd: LDT, a_win: NPFloats, l_win: NPFloats,
                 top_a: NPFloats, top_l: NPFloats, axes: DataAxes) -> Figure:
    """对真实采样做三次样条加密后调用 render。"""
    n: int = len(wd)
    sub: int = runtime.CFG["render"]["sub_per_day"]
    x: NPFloats = numpy.arange(n, dtype=float)
    xf: NPFloats = numpy.linspace(x[0], x[-1], n * sub)
    a_smooth: NPFloats = numpy.clip(make_interp_spline(x, a_win, k=3)(xf), 0.0, None)
    l_smooth: NPFloats = numpy.clip(make_interp_spline(x, l_win, k=3)(xf), 0.0, None)
    return render(wd, a_win, l_win, top_a, top_l, xf, a_smooth, l_smooth,
                  axes["left_ticks"], axes["right_ticks"], float(axes["ymax"]),
                  float(axes["right_top"]), runtime.CFG["app"]["x_label"])

def save_figure(fig: Figure, out_path: str) -> None:
    """保存图表；输出格式由扩展名决定。"""
    fig.savefig(out_path, facecolor="white")
    matplotlib.pyplot.close(fig)
    log("      wrote %s (%d bytes)" % (out_path, os.path.getsize(out_path)))

def metadata_path(out_path: str) -> str:
    """与图片同名、同目录，扩展名替换为 .json。"""
    return os.path.splitext(out_path)[0] + ".json"

def write_metadata(out_path: str, wd: LDT, last_recorded: str, src_r: str, pages: DataPages,
                   window: int, a_win: NPFloats, l_win: NPFloats, axes: DataAxes,
                   v_ana: NPFloats, v_let: NPFloats, top_a: NPFloats, top_l: NPFloats,
                   info_a: DataRow, info_l: DataRow, ana: Drug, let: Drug):
    """写出与图片同名的 JSON 元数据。"""
    app: DStr = runtime.CFG["app"]
    c: DataRow = runtime.CFG["colors"]
    n_int: int = runtime.CFG["render"]["n_intervals"]
    mean_h: float = float(axes["mean_h"])

    meta: TMeta = {
        "file": os.path.basename(out_path),
        "generated_by": os.path.basename(__file__),
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
                           float(numpy.abs(numpy.diff(v_ana + v_let)).mean()), 4),
                       "plotted_mean_abs_day_change": round(
                           float(numpy.abs(numpy.diff(top_a + top_l)).mean()), 4),
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
        "pk": {"anastrozole": {"t_half_h": runtime.CFG["pk"]["anastrozole"]["t_half_h"],
                               "ka": runtime.CFG["pk"]["anastrozole"]["ka"],
                               "V_per_F_L": round(ana.V or 0, 2),
                               "calibration": "稳态谷值 25.7 ng/mL @1 mg/日"},
               "letrozole": {"t_half_h": runtime.CFG["pk"]["letrozole"]["t_half_h"],
                             "ka": runtime.CFG["pk"]["letrozole"]["ka"],
                             "V_per_F_L": round(let.V or 0, 2)}},
    }
    meta_path: str = metadata_path(out_path)
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    log("      -> %s" % meta_path)
