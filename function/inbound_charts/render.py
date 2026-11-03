import matplotlib, matplotlib.pyplot, numpy
from datetime import datetime
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
from matplotlib.artist import Artist
from type_alias import LStr, LInt, SStr, LDT, DStr, NPFloats, DataRow
import runtime

def setup_fonts() -> str:
    """按配置设置等宽字体，返回实际使用的字体名。"""
    from matplotlib.font_manager import fontManager
    prefer: LStr = runtime.CFG["fonts"]["monospace"]
    have: SStr = {f.name for f in fontManager.ttflist}
    picked: LStr = [f for f in prefer if f in have]
    matplotlib.rcParams["font.family"] = "monospace"
    matplotlib.rcParams["font.monospace"] = picked + ["DejaVu Sans Mono"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    return picked[0] if picked else "DejaVu Sans Mono"

def render(wd: LDT, a_win: NPFloats, l_win: NPFloats, b_ana: NPFloats,
           b_let: NPFloats, xf: NPFloats, a_smooth: NPFloats, l_smooth: NPFloats,
           left_ticks: float|NPFloats, right_ticks: float|NPFloats,
           ymax: float, right_top: float, xlabel: str) -> Figure:
    """生成双轴图表。"""
    c: DataRow = runtime.CFG["colors"]
    app: DStr = runtime.CFG["app"]
    n_int: int = runtime.CFG["render"]["n_intervals"]

    fig: Figure; ax: Axes; ax2: Axes
    fig, ax = matplotlib.pyplot.subplots(figsize=(16.0, 9.0), dpi=110)
    fig.patch.set_facecolor("white")
    ax2 = ax.twinx()
    ax2.set_zorder(ax.get_zorder() - 1)
    ax.patch.set_visible(False)

    xs: NPFloats = numpy.linspace(xf[0], xf[-1], len(b_ana))
    ax2.fill_between(xs, 0.0, b_ana, color=c["ana_bar"], linewidth=0, zorder=2)
    ax2.fill_between(xs, b_ana, b_ana + b_let, color=c["let_bar"], linewidth=0, zorder=2)
    ax2.set_ylim(0, right_top)

    ax.plot(xf, a_smooth, color=c["ana_line"], linewidth=2.6,
            solid_capstyle="round", label=app["label_anastrozole"], zorder=5)
    ax.plot(xf, l_smooth, color=c["let_line"], linewidth=2.6,
            solid_capstyle="round", label=app["label_letrozole"], zorder=5)
    x: NPFloats = numpy.arange(len(wd), dtype=float)
    ax.scatter(x, a_win, s=7, color=c["ana_line"], alpha=0.28, zorder=4)
    ax.scatter(x, l_win, s=7, color=c["let_line"], alpha=0.28, zorder=4)

    ax.set_ylim(0, ymax)
    ax.set_xlim(-0.6, len(wd) - 0.4)
    ax.set_yticks(left_ticks)
    ax2.set_yticks(right_ticks)
    # 右轴标签小数位随步长自适应，避免小范围图表退化为「0 / 1」
    r_step: float = right_top / n_int
    r_dec: int = 0 if r_step >= 1.0 else (2 if r_step < 0.1 else 1)
    ticks_arr: NPFloats = numpy.atleast_1d(right_ticks)
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
