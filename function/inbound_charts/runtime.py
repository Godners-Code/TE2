import os
from type_alias import TConfig

CFG: TConfig
SCRIPT_BASENAME: str

def default_config_name() -> str:
    script_stem: str = os.path.splitext(SCRIPT_BASENAME)[0]
    return script_stem + ".toml"

# V_per_F == 0 时必须同时提供的标定参数
PK_CALIB_KEYS: tuple[str, str, str] = ("calib_trough", "calib_dose", "calib_tau")

def log(*content: str) -> None:
    """统一日志输出（立即刷新）。"""
    print(*content, flush=True)
