# -*- coding: utf-8 -*-
"""
T/E2 ADJUST — one-shot inbound pipeline
=======================================
Whole job in a single script:

  1. fetch README.html and resolve the newest "逐日明细" page, i.e. the
     Inbound-<year>.html link with the largest year (falls back to
     Inbound-<year-1>.html when the 6-month window crosses a year boundary);
  2. parse the daily records and take the anastrozole / letrozole doses of the
     last 6 months (180 days);
  3. run the one-compartment first-order-absorption PK model (linear
     superposition of every dose) and draw the finalised dual-axis chart.

Output: out\\Inbound-<end-date>.png  +  out\\chart-meta.json

Usage
-----
  python te2_chart_pipeline.py [--end auto|YYYY-MM-DD] [--days 180]
                               [--window 20] [--no-fetch] [--outdir DIR]

Notes
-----
* Network: direct schannel/curl TLS handshakes fail on this host, so the fetch
  goes through urllib + the local proxy; pages are cached under te2-site/.
* The 08:00 daily sample is the pre-dose trough, matching the settled spec.
"""
import argparse
import csv
import html as html_mod
import json
import math
import os
import re
import ssl
import sys
import time
import urllib.request
from datetime import datetime, timedelta

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline, make_interp_spline

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------
BASE = r"D:\Manuals\TE2"
SITE = "https://te2.tianyue.ren/"
CACHE = os.path.join(BASE, "te2-site")
OUT = os.path.join(BASE, "out")
PROXY = "http://192.168.109.5:8080"

SAMPLE_HOUR = 8            # dose AND sampling clock hour -> pre-dose trough
SAMPLE_STEP_H = 24.0
SUB_PER_DAY = 14
N_INT = 6                  # equal interval count on both y axes
DOSE_LOOKBACK_H = 24 * 40  # beyond ~20 half-lives a dose no longer matters
CONC_DECIMALS = 4          # canonical series precision, same as the CSV store

ANA_LINE, LET_LINE = "#1a73e8", "#12b76a"
ANA_BAR, LET_BAR = "#a9cbe8", "#a6d9b8"
GRID_COLOR, GRID_LW = "#b0bac4", 0.9

TITLE = "T/E2 ADJUST Inbound Deduce"
YLABEL_LEFT = "Blood Concentration (ng/mL)"
YLABEL_RIGHT = "Volatility on Last 20 Days (Standard Deviation)"
XLABEL = "Date (Daily Sampling)"
anastrozole_label = "Anastrozole"
letrozole_label = "Letrozole"

# calibrated PK parameters (see pk_te2.py for the calibration record)
PK = {
    "anastrozole": {"t_half_h": 50.0, "ka": 2.5, "V_per_F": None,   # calibrated below
                    "calib_trough": 25.7, "calib_dose": 1.0, "calib_tau": 24.0},
    "letrozole":   {"t_half_h": 48.0, "ka": 2.0, "V_per_F": 1.9 * 70.0},
}


def log(*a):
    print(*a, flush=True)


# --------------------------------------------------------------------------
# 1. fetching
# --------------------------------------------------------------------------
def _opener():
    ctx = ssl.create_default_context()
    return urllib.request.build_opener(
        urllib.request.ProxyHandler({"http": PROXY, "https": PROXY}),
        urllib.request.HTTPSHandler(context=ctx))


def fetch(name, use_network=True, timeout=60):
    """fetch <name> from the site, caching it under te2-site/.
    Falls back to the cached copy when the network is unavailable."""
    path = os.path.join(CACHE, name)
    if use_network:
        url = SITE + name
        t0 = time.time()
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with _opener().open(req, timeout=timeout) as r:
                raw = r.read()
            os.makedirs(CACHE, exist_ok=True)
            with open(path, "wb") as f:
                f.write(raw)
            log("      GET %-24s %6d bytes  %.2fs" % (name, len(raw), time.time() - t0))
            return raw.decode("utf-8", "replace"), "network"
        except Exception as e:
            log("      GET %-24s FAILED (%s: %s)" % (name, type(e).__name__, e))
    if os.path.exists(path):
        txt = open(path, encoding="utf-8", errors="replace").read()
        log("      CACHE %-22s %6d bytes" % (name, os.path.getsize(path)))
        return txt, "cache"
    return None, "missing"


def newest_inbound_year(readme_html):
    """resolve the newest Inbound-<year>.html link advertised by README"""
    years = sorted({int(y) for y in re.findall(r'Inbound-(\d{4})\.html', readme_html)})
    if not years:
        raise SystemExit("README does not link any Inbound-<year>.html page")
    return years[-1], years


# --------------------------------------------------------------------------
# 2. parsing the daily records
# --------------------------------------------------------------------------
def cell_text(td):
    td = re.sub(r"(?is)<br\s*/?>", " / ", td)
    td = re.sub(r"(?s)<[^>]+>", " ", td)
    td = html_mod.unescape(td).replace("\u00a0", " ")
    return re.sub(r"\s+", " ", td).strip()


def num(s):
    s = (s or "").strip()
    if s in ("", "-", "--"):
        return None
    return float(s) if re.match(r"^-?\d+(?:\.\d+)?$", s) else None


def parse_inbound(page_html):
    """-> list of {date, anastrozole_mg, letrozole_mg, ...} sorted by date"""
    rows, seen = [], set()
    for tr in re.findall(r"(?is)<tr>(.*?)</tr>", page_html):
        tds = re.findall(r"(?is)<td[^>]*>(.*?)</td>", tr)
        if len(tds) < 6:
            continue
        v = [cell_text(c) for c in tds]
        if not re.match(r"^\d{4}-\d{2}-\d{2}$", v[0]) or v[0] in seen:
            continue
        seen.add(v[0])
        rows.append({"date": v[0],
                     "anastrozole_mg": num(v[2]),
                     "letrozole_mg": num(v[3]),
                     "testosterone_mg": num(v[4]) if len(v) > 4 else None,
                     "hcg_IU": num(v[5]) if len(v) > 5 else None})
    rows.sort(key=lambda r: r["date"])
    return rows


# --------------------------------------------------------------------------
# 3. PK model (one compartment, first-order absorption, linear superposition)
# --------------------------------------------------------------------------
class Drug:
    def __init__(self, key):
        p = dict(PK[key])
        self.name = key
        self.ke = math.log(2.0) / p["t_half_h"]
        self.ka = p["ka"]
        self.V = p["V_per_F"]
        if self.V is None:                      # calibrate to the label trough
            R = math.exp(-self.ke * p["calib_tau"])
            shape = (self.ka / (self.ka - self.ke)) * \
                    (math.exp(-self.ke * p["calib_tau"]) - math.exp(-self.ka * p["calib_tau"]))
            self.V = 1000.0 * p["calib_dose"] * shape / (1.0 - R) / p["calib_trough"]

    def conc(self, dt_h, dose_mg):
        if dose_mg is None or dose_mg <= 0 or dt_h < 0:
            return 0.0
        c = (dose_mg * self.ka) / (self.V * (self.ka - self.ke))
        return 1000.0 * c * (math.exp(-self.ke * dt_h) - math.exp(-self.ka * dt_h))

    def concentration(self, when, doses):
        """sum the contributions of every dose still within the lookback window.

        The lookback cutoff matters: a 1 mg dose older than ~20 half-lives still
        contributes ~1e-5 ng/mL, which is enough to flip the 4th decimal (and
        therefore an anti-aliased pixel) relative to a simulation that truncates.
        """
        total = 0.0
        for t, mg in doses:
            dt_h = (when - t).total_seconds() / 3600.0
            if dt_h < 0 or dt_h > DOSE_LOOKBACK_H:
                continue
            total += self.conc(dt_h, mg)
        return total


def build_doses(rows, field):
    return [(datetime.strptime(r["date"], "%Y-%m-%d") + timedelta(hours=SAMPLE_HOUR), r[field])
            for r in rows if r.get(field)]


def daily_series(drug, doses, start, end):
    """one sample per day at SAMPLE_HOUR, inclusive of both ends"""
    out, t = [], start
    while t <= end:
        out.append(drug.concentration(t, doses))
        t += timedelta(hours=SAMPLE_STEP_H)
    return np.array(out)


def canon(series, decimals=CONC_DECIMALS):
    """round exactly the way the CSV store does: format to a decimal string and
    parse it back, so the resulting doubles are the ones a CSV round-trip yields
    (np.round can land 1 ULP away and shift anti-aliased pixels)."""
    return np.array([float(("%." + str(decimals) + "f") % v) for v in series])


def rolling_vol(x, window):
    out = np.full(len(x), np.nan)
    for i in range(window - 1, len(x)):
        out[i] = x[i - window + 1: i + 1].std(ddof=1)
    return out


def fit_bar_curve(x, y):
    sigma = float(np.std(np.diff(y), ddof=1) / np.sqrt(2.0))
    s = len(x) * sigma ** 2
    spl = UnivariateSpline(x, y, k=3, s=s)
    return spl, {"sigma_est": round(sigma, 6), "s": round(s, 6),
                 "residual_rms": round(float(np.sqrt(np.mean((y - np.clip(spl(x), 0, None)) ** 2))), 6),
                 "raw_mean": round(float(y.mean()), 6),
                 "fitted_mean": round(float(np.clip(spl(x), 0, None).mean()), 6)}


# --------------------------------------------------------------------------
# 4. rendering
# --------------------------------------------------------------------------
def setup_fonts():
    from matplotlib import font_manager
    prefer = ["Microsoft YaHei", "Microsoft YaHei UI", "SimHei", "SimSun"]
    have = {f.name for f in font_manager.fontManager.ttflist}
    picked = [f for f in prefer if f in have]
    if picked:
        matplotlib.rcParams["font.sans-serif"] = picked + ["DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    return picked[0] if picked else "NONE"


def render(wd, a_win, l_win, b_ana, b_let, xf, a_smooth, l_smooth,
           left_ticks, right_ticks, ymax, right_top, xlabel):
    fig, ax = plt.subplots(figsize=(16.0, 9.0), dpi=110)
    fig.patch.set_facecolor("white")

    ax2 = ax.twinx()
    ax2.set_zorder(ax.get_zorder() - 1)
    ax.patch.set_visible(False)

    xs = np.linspace(xf[0], xf[-1], len(b_ana))
    ax2.fill_between(xs, 0.0, b_ana, color=ANA_BAR, linewidth=0, zorder=2)
    ax2.fill_between(xs, b_ana, b_ana + b_let, color=LET_BAR, linewidth=0, zorder=2)
    ax2.set_ylim(0, right_top)

    ax.plot(xf, a_smooth, color=ANA_LINE, linewidth=2.6, solid_capstyle="round",
            label=anastrozole_label, zorder=5)
    ax.plot(xf, l_smooth, color=LET_LINE, linewidth=2.6, solid_capstyle="round",
            label=letrozole_label, zorder=5)
    x = np.arange(len(wd), dtype=float)
    ax.scatter(x, a_win, s=7, color=ANA_LINE, alpha=0.28, zorder=4)
    ax.scatter(x, l_win, s=7, color=LET_LINE, alpha=0.28, zorder=4)

    ax.set_ylim(0, ymax)
    ax.set_xlim(-0.6, len(wd) - 0.4)
    ax.set_yticks(left_ticks)
    ax2.set_yticks(right_ticks)
    # label precision follows the right-axis step, so a small-range chart does not
    # collapse to "0 / 1" (the default step ~5.0 still yields plain integers)
    r_step = right_top / N_INT
    r_dec = 0 if r_step >= 1.0 else (2 if r_step < 0.1 else 1)
    ax2.set_yticklabels([("%." + str(r_dec) + "f") % v for v in right_ticks])

    ax.set_ylabel(YLABEL_LEFT, fontsize=12.5, color="#1a202c")
    ax2.set_ylabel(YLABEL_RIGHT, fontsize=12.5, color="#4a5568")
    ax.tick_params(axis="y", colors="#1a202c", labelsize=10)
    ax2.tick_params(axis="y", colors="#4a5568", labelsize=10)

    tick_idx = [i for i, d in enumerate(wd) if d.day in (1, 16)]
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([wd[i].strftime("%m-%d") for i in tick_idx], fontsize=10)
    ax.set_xlabel(xlabel, fontsize=11.5)

    for i in tick_idx:
        ax.axvline(i, color=GRID_COLOR, linewidth=GRID_LW, zorder=3)
    for v in left_ticks[1:]:
        ax.axhline(v, color=GRID_COLOR, linewidth=GRID_LW, zorder=3)
    ax.set_axisbelow(False)

    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax.spines["left"].set_color("#a0aec0")
    ax2.spines["right"].set_color("#a0aec0")

    ax.set_title(TITLE, fontsize=16, fontweight="bold", pad=14)
    h1, l1 = ax.get_legend_handles_labels()
    ax.legend(h1, l1, loc="upper left", fontsize=10.5, ncol=2, framealpha=0.0)

    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--end", default="auto", help="'auto' = last recorded day")
    ap.add_argument("--days", type=int, default=180, help="≈ last 6 months")
    ap.add_argument("--window", type=int, default=20)
    ap.add_argument("--no-fetch", action="store_true")
    ap.add_argument("--outdir", default=OUT)
    args = ap.parse_args()
    use_net = not args.no_fetch
    os.makedirs(args.outdir, exist_ok=True)

    log("[1/6] resolving the newest daily-record page from README.html ...")
    readme, src_r = fetch("README.html", use_net)
    if readme is None:
        raise SystemExit("README.html unavailable (network and cache both failed)")
    year, all_years = newest_inbound_year(readme)
    log("      README source=%s  Inbound years seen=%s  -> newest = %d"
        % (src_r, all_years, year))

    log("[2/6] fetching and parsing daily records ...")
    pages = {}
    target = "Inbound-%d.html" % year
    txt, src = fetch(target, use_net)
    if txt is None:
        raise SystemExit("%s unavailable" % target)
    pages[year] = parse_inbound(txt)
    log("      %s (%s): %d daily rows  %s .. %s"
        % (target, src, len(pages[year]), pages[year][0]["date"], pages[year][-1]["date"]))

    rows = list(pages[year])
    last_recorded = next(r["date"] for r in reversed(rows)
                         if r["anastrozole_mg"] or r["letrozole_mg"])
    end_date = (datetime.strptime(last_recorded, "%Y-%m-%d").date()
                if args.end == "auto" else datetime.strptime(args.end, "%Y-%m-%d").date())
    log("      last recorded dosing day = %s  -> window end = %s" % (last_recorded, end_date))

    log("[3/6] extending the window backwards ...")
    win_end = datetime.combine(end_date, datetime.min.time()) + timedelta(hours=SAMPLE_HOUR)
    win_start = win_end - timedelta(days=args.days - 1)
    sim_start = win_start - timedelta(days=args.window - 1 + 40)   # rolling window + PK lookback
    log("      window    %s .. %s  (%d days)"
        % (win_start.date(), win_end.date(), args.days))
    log("      simulated %s .. %s  (extra %d days for the rolling window and dose lookback)"
        % (sim_start.date(), win_end.date(), (win_start - sim_start).days))

    if sim_start.year < year:
        prev = "Inbound-%d.html" % (year - 1)
        t2, s2 = fetch(prev, use_net)
        if t2 is None:
            log("      %s not available - doses before %d-01-01 are unknown; "
                "the first weeks may under-read" % (prev, year))
        else:
            pages[year - 1] = parse_inbound(t2)
            rows = pages[year - 1] + rows
            log("      %s (%s): %d daily rows %s .. %s"
                % (prev, s2, len(pages[year - 1]),
                   pages[year - 1][0]["date"], pages[year - 1][-1]["date"]))
            rows.sort(key=lambda r: r["date"])

    doses = {k: build_doses(rows, k) for k in ("anastrozole_mg", "letrozole_mg")}
    for k, v in doses.items():
        log("      %-16s doses=%3d  total=%7.1f mg  first=%s"
            % (k, len(v), sum(m for _, m in v), v[0][0].date() if v else "-"))

    log("[4/6] PK simulation (1 point/day @%02d:00, pre-dose trough) ..." % SAMPLE_HOUR)
    ana, let = Drug("anastrozole"), Drug("letrozole")
    log("      anastrozole t1/2=%.1fh ka=%.2f V/F=%.1f L | letrozole t1/2=%.1fh ka=%.2f V/F=%.1f L"
        % (ana.ke and math.log(2) / ana.ke, ana.ka, ana.V,
           math.log(2) / let.ke, let.ka, let.V))
    a_full = canon(daily_series(ana, doses["anastrozole_mg"], sim_start, win_end))
    l_full = canon(daily_series(let, doses["letrozole_mg"], sim_start, win_end))
    log("      simulated points=%d" % len(a_full))

    off = (win_start - sim_start).days
    a_win, l_win = a_full[off:], l_full[off:]
    n = len(a_win)
    x = np.arange(n, dtype=float)
    log("      window points=%d  peak ana=%.2f ng/mL  peak let=%.2f ng/mL"
        % (n, a_win.max(), l_win.max()))

    log("[5/6] rolling %d-day volatility + smooth band top ..." % args.window)
    v_ana = np.nan_to_num(rolling_vol(a_full, args.window))[off:]
    v_let = np.nan_to_num(rolling_vol(l_full, args.window))[off:]
    spl_a, info_a = fit_bar_curve(x, v_ana)
    spl_l, info_l = fit_bar_curve(x, v_let)
    xs_sub = np.linspace(x[0], x[-1], n * SUB_PER_DAY)
    top_a = np.clip(spl_a(xs_sub), 0.0, None)
    top_l = np.clip(spl_l(xs_sub), 0.0, None)

    mean_h = float((top_a + top_l).mean())
    right_top = mean_h / 0.25
    data_max = float(max(a_win.max(), l_win.max()))
    step_l = float(math.ceil(data_max / N_INT))
    ymax = step_l * N_INT
    left_ticks = np.linspace(0.0, ymax, N_INT + 1)
    right_ticks = np.linspace(0.0, right_top, N_INT + 1)
    log("      fitted bar top |day-to-day|: raw=%.4f -> plotted=%.4f"
        % (np.abs(np.diff(v_ana + v_let)).mean(), np.abs(np.diff(top_a + top_l)).mean()))
    log("      mean=%.4f -> right axis 0..%.4f (%.1f%%)  left axis 0..%.0f step=%.0f"
        % (mean_h, right_top, 100 * mean_h / right_top, ymax, step_l))

    log("[6/6] rendering ...")
    font = setup_fonts()
    log("      CJK font -> %s" % font)
    xf = np.linspace(x[0], x[-1], n * SUB_PER_DAY)
    a_smooth = np.clip(make_interp_spline(x, a_win, k=3)(xf), 0.0, None)
    l_smooth = np.clip(make_interp_spline(x, l_win, k=3)(xf), 0.0, None)
    xlabel = XLABEL
    fig = render([win_start + timedelta(days=i) for i in range(n)], a_win, l_win,
                 top_a, top_l, xf, a_smooth, l_smooth,
                 left_ticks, right_ticks, ymax, right_top, xlabel)

    wd = [win_start + timedelta(days=i) for i in range(n)]
    name = "Inbound-%s.png" % wd[-1].strftime("%Y-%m-%d")
    path = os.path.join(args.outdir, name)
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    log("      wrote %s (%d bytes)" % (path, os.path.getsize(path)))

    meta = {
        "file": name,
        "generated_by": "tools/te2_chart_pipeline.py",
        "source": {"readme": src_r, "inbound_pages": sorted(pages.keys()),
                   "window_end_rule": "'auto' = last recorded dosing day (%s)" % last_recorded},
        "window": {"start": wd[0].strftime("%Y-%m-%d"), "end": wd[-1].strftime("%Y-%m-%d"),
                   "days": n, "sampling": "1 point per day @08:00 (pre-dose trough)"},
        "title": TITLE, "x_label": xlabel,
        "left_axis": {"label": YLABEL_LEFT, "max": round(ymax, 3),
                      "line": "interpolating cubic spline through every real daily sample",
                      "colors": {"anastrozole": ANA_LINE, "letrozole": LET_LINE},
                      "peaks": {"anastrozole": round(float(a_win.max()), 3),
                                "letrozole": round(float(l_win.max()), 3)}},
        "right_axis": {"label": YLABEL_RIGHT, "max": round(right_top, 4),
                       "mean_at_fraction": round(mean_h / right_top, 4),
                       "band_colors": {"anastrozole": ANA_BAR, "letrozole": LET_BAR},
                       "bar_top": "smooth fitted curve replacing the horizontal bar-top segments",
                       "moving_average_applied": False, "legend": False},
        "volatility": {"window_days": args.window, "statistic": "标准差",
                       "statistic_key": "std", "unit": "ng/mL",
                       "raw_mean_total": round(float((v_ana + v_let).mean()), 4),
                       "plotted_mean_total": round(mean_h, 4),
                       "raw_mean_abs_day_change": round(float(np.abs(np.diff(v_ana + v_let)).mean()), 4),
                       "plotted_mean_abs_day_change": round(float(np.abs(np.diff(top_a + top_l)).mean()), 4),
                       "fit_anastrozole": info_a, "fit_letrozole": info_l},
        "grid": {"vertical": "1st and 16th of every month",
                 "horizontal": "both y axes, identical interval count",
                 "color": GRID_COLOR, "linewidth": GRID_LW, "zorder": 3},
        "axis_intervals": {"intervals_per_axis": N_INT,
                           "left_step": round(step_l, 4), "left_max": round(ymax, 4),
                           "right_step": round(right_top / N_INT, 4),
                           "right_max": round(right_top, 4)},
        "mean_reference_line": False,
        "stage_extrema": {"drawn": False, "note": "no extrema markers in the finalised chart"},
        "pk": {"anastrozole": {"t_half_h": 50.0, "ka": 2.5, "V_per_F_L": round(ana.V, 2),
                               "calibration": "steady-state trough 25.7 ng/mL @1 mg/day"},
               "letrozole": {"t_half_h": 48.0, "ka": 2.0, "V_per_F_L": round(let.V, 2)}},
    }
    with open(os.path.join(args.outdir, "chart-meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    log("      -> chart-meta.json")
    log("DONE")


if __name__ == "__main__":
    main()
