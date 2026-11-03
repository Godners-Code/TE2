import os, tomllib, re
from type_alias import DFloat, LStr, DInt, DStr, DictParam, TConfig
from runtime import default_config_name, PK_CALIB_KEYS
import runtime

def default_config_path() -> str:
    """默认配置文件：与本脚本同目录的 <脚本名去掉扩展名>.toml。"""
    here: str = os.path.dirname(runtime.SCRIPT_BASENAME)
    return os.path.join(here, default_config_name())

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
    p: DInt = runtime.CFG["period"]
    if s in ("monthly", "quarter", "half", "year"): return int(p[s])
    if re.fullmatch(r"\d+", s): return int(s)
    raise SystemExit("配置错误：[run].period 必须为 monthly、quarter、half、"
                     "year 或整数天数，实际为 %r" % value)

def load_run_params() -> DictParam:
    """从 [run] / [output] / [site] 读取并校验全部运行参数。"""
    run:  DictParam = runtime.CFG["run"]
    out:  DStr      = runtime.CFG["output"]
    site: DStr      = runtime.CFG["site"]
    p:    DInt      = runtime.CFG["period"]

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
