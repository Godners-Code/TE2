import math, numpy
from datetime import datetime, timedelta
from scipy.interpolate import UnivariateSpline
from type_alias import DFloat, DataDoses, DataRows, DataRow, NPFloats, LFloat
import runtime

class Drug:
    """一室一级吸收模型。"""
    def __init__(self, key: str) -> None:
        p: DFloat = dict(runtime.CFG["pk"][key])
        self.name: str = key
        self.ke: float = math.log(2.0) / p["t_half_h"]
        self.ka: float = p["ka"]
        self.V: float|None = float(p.get("V_per_F", 0.0) or 0.0) or None
        if self.V is None:
            # V_per_F == 0：按说明书谷值反算 V/F
            # （标定参数的完整性在 load_config 中已校验）
            r: float = math.exp(-self.ke * p["calib_tau"])
            shape: float = (self.ka / (self.ka - self.ke)) * \
                    (math.exp(-self.ke * p["calib_tau"]) -
                     math.exp(-self.ka * p["calib_tau"]))
            self.V = 1000.0 * p["calib_dose"] * shape / (1.0 - r) / p["calib_trough"]

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
        lookback: int = runtime.CFG["render"]["dose_lookback_h"]
        for t, mg in doses:
            dt_h: float = (when - t).total_seconds() / 3600.0
            if dt_h < 0 or dt_h > lookback:
                continue
            total += self.conc(dt_h, float(mg))
        return total

def build_doses(rows: DataRows, field: str) -> DataDoses:
    """从每日记录中提取给药时间点与剂量。"""
    hour: int = runtime.CFG["render"]["sample_hour"]
    return [(datetime.strptime(str(r["date"]), "%Y-%m-%d") + timedelta(hours=hour), r[field])
            for r in rows if r.get(field)]

def daily_series(drug: Drug, doses: DataDoses,
                 start: datetime, end: datetime) -> NPFloats:
    """按每天一次、在 sample_hour 采样，首尾都包含。"""
    step: float = runtime.CFG["render"]["sample_step_h"]
    out: LFloat = []
    t: datetime = start
    while t <= end:
        out.append(drug.concentration(t, doses))
        t += timedelta(hours=step)
    return numpy.array(out)

def canon(series: NPFloats, decimals: int|None = None) -> NPFloats:
    """按 CSV 存储的方式保留小数：格式化成字符串再解析回来，
    保证结果与 CSV 往返后完全一致（numpy.round 可能偏移 1 ULP 并影响抗锯齿像素）。"""
    if decimals is None: decimals = runtime.CFG["render"]["conc_decimals"]
    fmt: str = "%." + str(decimals) + "f"
    return numpy.array([float(fmt % v) for v in series])

def rolling_vol(x: NPFloats, window: int) -> NPFloats:
    """滚动窗口标准差。"""
    out: NPFloats = numpy.full(len(x), numpy.nan)
    for i in range(window - 1, len(x)):
        out[i] = x[i - window + 1: i + 1].std(ddof=1)
    return out

def fit_bar_curve(x: NPFloats, y: NPFloats) -> tuple[UnivariateSpline, DataRow]:
    """按差分的经验标准差选取平滑参数，拟合色带顶部曲线。"""
    sigma: float = numpy.std(numpy.diff(y), ddof=1) / numpy.sqrt(2.0)
    s: float = len(x) * sigma ** 2
    spl: UnivariateSpline = UnivariateSpline(x, y, k=3, s=s)
    return spl, {"sigma_est": round(sigma, 6), "s": round(s, 6),
                 "residual_rms": round(numpy.sqrt(
                     numpy.mean((y - numpy.clip(spl(x), 0, None)) ** 2)), 6),
                 "raw_mean": round(float(y.mean()), 6),
                 "fitted_mean": round(float(numpy.clip(spl(x), 0, None).mean()), 6)}
