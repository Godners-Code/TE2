import argparse
from argparse import ArgumentParser, Namespace
from load_config import default_config_path
from runtime import default_config_name

def parse_args() -> Namespace:
    """命令行只保留 --config；其余参数全部来自 TOML。"""
    ap: ArgumentParser = argparse.ArgumentParser(
        description="TE2 入站图表流水线（所有参数见 TOML 配置文件）")
    ap.add_argument("--config", default=default_config_path(),
                    help="TOML 配置文件路径；默认与本脚本同目录的 %s"
                         % default_config_name())
    return ap.parse_args()
