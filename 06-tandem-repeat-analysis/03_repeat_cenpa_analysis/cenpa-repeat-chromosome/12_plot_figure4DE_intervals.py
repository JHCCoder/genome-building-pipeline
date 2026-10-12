#!/usr/bin/env python3
"""
Figure 4D/E, step 2: one small-multiples block over the chromosome-aligned intervals
around the dominant CENP-A-enriched region (replaces the former panels D, E and F).

  row D: CENP-A and H3K27ac mean coverage (1/NH-weighted), HiFi depth (all reads,
         primary alignments only)
  row E: % of interval bases covered by merged arrays of four period bins
         (1-10, 193-195, 348-349, 386-390 bp; Figure 4 family colours)

All eight plots share one x-axis (labelled once, on the bottom row); each has its own
y-axis. Shorter- and longer-arm flanks are shaded, with one key. Drawn at printed size
for Genome Research (8 pt text, no titles).

Reads $FIG4_INTERVALS_TSV from 11_prepare_figure4DE_intervals.py.

Run with the python-visualizations env:
    python 12_plot_figure4DE_intervals.py [--width-mm 175] [--height-mm 80]
"""

import argparse
import os
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
from matplotlib.patches import Patch


def _cfg(name):
    val = os.environ.get(name)
    if val:
        return val
    root = os.path.abspath(os.path.dirname(__file__))
    # repo-root config.sh (next to .git): this tree also has a local config.sh
    while not (os.path.isfile(os.path.join(root, "config.sh")) and os.path.exists(os.path.join(root, ".git"))) \
            and root != os.path.dirname(root):
        root = os.path.dirname(root)
    # Source config.sh in a subshell so nested $VAR references are expanded.
    val = subprocess.check_output(
        ["bash", "-c", f'source "$1" 2>/dev/null; printf "%s" "${{{name}}}"', "_",
         os.path.join(root, "config.sh")],
        text=True,
    )
    if not val:
        sys.exit(f"ERROR: {name} is not set in the environment or in config.sh")
    return val


FONT_SIZE = 8                      # pt at final print size; Genome Research: 8-10 pt
SHORT_TINT, LONG_TINT = "#E3EEF6", "#FBF1DC"
ROWS = [  # (row y-axis title, [(metric, label, colour, tick formatter)])
    ("Mean coverage", [
        ("CENP-A", "CENP-A", "#1B4F72", "k"),
        ("H3K27ac", "H3K27ac", "#E69F00", "k"),
        ("HiFi all reads", "HiFi, all reads (×)", "#000000", None),
        ("HiFi primary", "HiFi, primary (×)", "#000000", None),
    ]),
    ("Bases in arrays (%)", [
        ("1-10 bp", "1–10 bp", _cfg("FIG4_COLOR_1_10"), None),
        ("193-195 bp", "193–195 bp", _cfg("FIG4_COLOR_193_195"), None),
        ("348-349 bp", "348–349 bp", _cfg("FIG4_COLOR_348_349"), None),
        ("386-390 bp", "386–390 bp", _cfg("FIG4_COLOR_386_390"), None),
    ]),
]
XLABELS = ["2\u20135", "0.5\u20132", "<0.5", "Core", "<0.5", "0.5\u20132", "2\u20135", "Background"]   # Mb, in the axis title

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=175)
ap.add_argument("--height-mm", type=float, default=80)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4

d = pd.read_csv(_cfg("FIG4_INTERVALS_TSV"), sep="\t")
out_dir = _cfg("FIG4_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, "figure4DE_interval_profiles")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 0.50, 0.05, 0.22, 0.74
COL_GAP_IN, ROW_GAP_IN = 0.36, 0.14
ncol = 4
cw = (PANEL_W - LEFT_IN - RIGHT_IN - (ncol - 1) * COL_GAP_IN) / ncol
rh = (PANEL_H - TOP_IN - BOTTOM_IN - ROW_GAP_IN) / 2
fig = plt.figure(figsize=(PANEL_W, PANEL_H))


def fmt_k(v, _):
    return "0" if v == 0 else (f"{v/1000:g}k" if abs(v) >= 1000 else f"{v:g}")


axes_rows = []
row_titles = []   # placed after drawing, 1 mm left of each row's widest y tick label
for r, (ytitle, plots) in enumerate(ROWS):
    y0 = PANEL_H - TOP_IN - (r + 1) * rh - r * ROW_GAP_IN
    row_axes = []
    for c, (metric, label, color, fmt) in enumerate(plots):
        ax = fig.add_axes(((LEFT_IN + c * (cw + COL_GAP_IN)) / PANEL_W, y0 / PANEL_H, cw / PANEL_W, rh / PANEL_H))
        s = d[d["metric"] == metric].sort_values("x")
        ax.axvspan(-0.5, 2.5, color=SHORT_TINT, lw=0, zorder=0)
        ax.axvspan(3.5, 6.5, color=LONG_TINT, lw=0, zorder=0)
        ax.plot(s["x"], s["value"], color=color, lw=1.2, marker="o", ms=2.6, zorder=3)
        ax.set_xlim(-0.5, 7.5)
        ax.set_ylim(0, s["value"].max() * 1.22)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(3, min_n_ticks=3))
        if fmt == "k":
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(fmt_k))
        ax.set_xticks(range(8))
        if r == len(ROWS) - 1:
            ax.set_xticklabels(XLABELS, rotation=45, ha="right", rotation_mode="anchor")
        else:
            ax.set_xticklabels([])
        ax.text(0.03, 0.97, label, transform=ax.transAxes, ha="left", va="top",
                bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.85), zorder=4)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        row_axes.append(ax)
    row_titles.append((fig.text(0, (y0 + rh / 2) / PANEL_H, ytitle, rotation=90, ha="right", va="center"),
                       row_axes[0]))
    axes_rows.append(row_axes)

fig.text((LEFT_IN + (PANEL_W - LEFT_IN - RIGHT_IN) / 2) / PANEL_W, 0.03 / PANEL_H,
         "Distance from the dominant CENP-A–enriched region (Mb)", ha="center", va="bottom")
fig.legend(handles=[Patch(color=SHORT_TINT, label="Shorter arm"), Patch(color=LONG_TINT, label="Longer arm")],
           ncol=2, loc="upper right", bbox_to_anchor=(1 - RIGHT_IN / PANEL_W, 1), frameon=False,
           handlelength=1.2, handleheight=0.9, columnspacing=1.0, borderaxespad=0.05)

fig.canvas.draw()
rend = fig.canvas.get_renderer()
for title, ax0 in row_titles:
    tick_left = min(t.get_window_extent(rend).x0 for t in ax0.get_yticklabels() if t.get_text()) / fig.dpi
    title.set_x((tick_left - 0.04) / PANEL_W)
    assert title.get_window_extent(rend).x0 >= 0, "row title runs off the left edge"
for ax in axes_rows[-1]:
    lbl_bottom = min(t.get_window_extent(rend).y0 for t in ax.get_xticklabels()) / fig.dpi
    assert lbl_bottom > 0.16, "x tick labels run into the axis title"

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
