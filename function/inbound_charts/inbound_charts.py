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
  python inbound_charts.py [--config PATH]

所有可调参数均来自 TOML 配置文件：
  [run].end / [run].period / [run].window
  [output].dir / [output].path
  [site].proxy
默认配置文件名为“脚本名去掉扩展名 + .toml”，
例如脚本为 inbound_charts.py 时，配置文件名是 inbound_charts.toml。
命令行只保留 --config，用于在需要时显式指定配置文件路径。

说明
----
* 需要 Python 3.12 或更高版本。
* 网络：通过 urllib 抓取；是否使用代理由 [site].proxy 决定
  （为空时显式禁用代理，非空时作为 HTTP/HTTPS 代理使用）。
* 页面缓存在系统临时目录：<temp>/<site.cache_dir_name>/。
* 每日 08:00 的采样为给药前谷值，与既定规范一致。
"""
# region # 引用库
import sys
if sys.version_info < (3, 12):
    raise SystemExit("Python 3.12 或更高版本是必需的。")
import os, numpy
from argparse import Namespace
from matplotlib.figure import Figure
from datetime import datetime, date, timedelta
from type_alias import LDT, NPFloats, DictParam, DataRow, DataRows, DataPages, DBDoses, DataAxes
from runtime import log
import runtime
runtime.SCRIPT_BASENAME = os.path.abspath(__file__)
from load_config import load_config, load_run_params
from pk_model import Drug
from render import setup_fonts
from output_path import resolve_output_path
from pipeline import (resolve_page_year, fetch_main_page, find_last_recorded,
                      resolve_end_date, compute_window, maybe_extend_previous_year,
                      build_dose_table, compute_bands, compute_axes, render_chart,
                      save_figure, write_metadata, simulate_pk)
from entry import parse_args
# endregion

# type_alias.py     类型别名
# runtime.py        运行时配置
# load_config.py    0. 配置加载
# fetch.py          1. 抓取
# parse_records.py  2. 解析记录
# pk_model.py       3. PK 模型
# render.py         4. 渲染
# output_path.py    5. 输出路径
# pipeline.py       6. 流水线各阶段
# entry.py          7. 解析参数

def inbound_charts_main():
    args: Namespace = parse_args()
    runtime.CFG = load_config(args.config)
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
        % runtime.CFG["render"]["sample_hour"])
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
        % (numpy.abs(numpy.diff(v_ana + v_let)).mean(),
           numpy.abs(numpy.diff(top_a + top_l)).mean()))

    axes: DataAxes = compute_axes(a_win, l_win, top_a, top_l)
    log("      mean=%.4f -> right axis 0..%.4f (%.1f%%)  left axis 0..%.0f step=%.0f"
        % (float(axes["mean_h"]), float(axes["right_top"]),
           100 * axes["mean_h"] / axes["right_top"],
           float(axes["ymax"]), float(axes["step_l"])))

    # 阶段 5：渲染与输出
    log("[6/6] rendering ...")
    font: str = setup_fonts()
    log("      monospace font -> %s" % font)

    wd: LDT = [win_start + timedelta(days=i) for i in range(n)]
    fig: Figure = render_chart(wd, a_win, l_win, top_a, top_l, axes)

    default_ext: str = runtime.CFG["image"]["default_ext"]
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

if __name__ == "__main__":
    inbound_charts_main()

