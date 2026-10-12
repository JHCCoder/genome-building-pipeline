#!/usr/bin/env python3
"""
Figure 4A — every TRF tandem repeat in the degu assembly: repeat-unit (monomer)
length versus total array length.

Same data and filtering as the Figure 4 version in check_tandem_repeats_degu.ipynb
(figure/tandem-repeat-investigation): chromosomes chr1-chr28 and chrX (unplaced
"seq" scaffolds and chrY dropped); point area scales with array length.

Revised for Genome Research (8 pt text, no title): all repeats grey; the three
CENP-A-relevant period families (193-195, 348-349, 386-390 bp) coloured and
labelled at the right instead of dashed lines plus a 29-entry chromosome key;
plain x tick labels (log10 superscripts would print below 8 pt). The point layer
is rasterized at 600 dpi (~1 million points); text and axes stay vector.

Run with the python-visualizations env:
    python 01_plot_figure4A_trf_scatter.py [--width-mm 56] [--height-mm 85]
"""

import argparse
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker


def _cfg(name):
    val = os.environ.get(name)
    if val:
        return val
    root = os.path.abspath(os.path.dirname(__file__))
    while not os.path.isfile(os.path.join(root, "config.sh")) and root != os.path.dirname(root):
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


FONT_SIZE = 8                     # pt at final print size; Genome Research: 8-10 pt
CHROMS = [f"chr{i}" for i in range(1, 29)] + ["chrX"]
FAMILIES = [  # (label, lowest period, highest period, colour); the axis already says bp
    ("193–195", 193, 195, _cfg("FIG4_COLOR_193_195")),
    ("348–349", 348, 349, _cfg("FIG4_COLOR_348_349")),
    ("386–390", 386, 390, _cfg("FIG4_COLOR_386_390")),
]
GREY = "#B0B0B0"
SIZE_MIN, SIZE_MAX = 0.3, 45      # pt^2, point area range (scaled by array length)

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=56)
ap.add_argument("--height-mm", type=float, default=85)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4

out_dir = _cfg("FIG4_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, "figure4A_trf_scatter")

df = pd.read_csv(_cfg("FIG4_TRF_TABLE"), sep="\t", usecols=["sequence", "start", "end", "period_size"])
df = df[df["sequence"].isin(CHROMS)].copy()
df["length"] = df["end"] - df["start"]
df = df[df["length"] > 0].sample(frac=1, random_state=28)       # shuffle so no chromosome sits on top
LMAX = df["length"].max()
area = lambda L: SIZE_MIN + (SIZE_MAX - SIZE_MIN) * L / LMAX
print(f"repeats: {len(df):,} on {df['sequence'].nunique()} chromosomes; longest array {LMAX:,} bp")
for lab, lo, hi, _ in FAMILIES:
    sub = df[df["period_size"].between(lo, hi)]
    print(f"  {lab}: {len(sub):,} repeats, {sub['length'].sum()/1e6:.1f} Mb")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "xtick.minor.width": 0.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 0.40, 0.58, 0.30, 0.40   # bottom matches 4B so the x-axes align
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
ax = fig.add_axes((LEFT_IN / PANEL_W, BOTTOM_IN / PANEL_H,
                   1 - (LEFT_IN + RIGHT_IN) / PANEL_W, 1 - (TOP_IN + BOTTOM_IN) / PANEL_H))

in_fam = np.zeros(len(df), bool)
for _, lo, hi, _ in FAMILIES:
    in_fam |= df["period_size"].between(lo, hi).values
bg = df[~in_fam]
ax.scatter(bg["length"], bg["period_size"], s=area(bg["length"]), c=GREY, lw=0, alpha=0.6, rasterized=True)
for lab, lo, hi, col in FAMILIES:
    sub = df[df["period_size"].between(lo, hi)]
    ax.scatter(sub["length"], sub["period_size"], s=area(sub["length"]), c=col, lw=0, alpha=0.85,
               rasterized=True, zorder=3)
    ax.text(1.03, (lo + hi) / 2, lab, transform=ax.get_yaxis_transform(), color=col,
            ha="left", va="center")

ax.set_xscale("log")
ax.set_xlim(df["length"].min() * 0.8, 2e7)
ax.set_xticks([1e2, 1e4, 1e6], ["100", "10k", "1M"])              # every other decade labelled
ax.set_xticks([1e3, 1e5, 1e7], minor=True)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylim(0, 520)
ax.set_yticks(range(0, 501, 100))
ax.set_xlabel("Total array length (bp)", labelpad=2)
ax.set_ylabel("Repeat-unit length (bp)", labelpad=2)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

# Size key above the plot, one line: "Array length  * 1 Mb  * 3 Mb  * 6 Mb"
key_ax = fig.add_axes((0, 1 - 0.24 / PANEL_H, 1, 0.20 / PANEL_H))
key_ax.set_axis_off()
key_ax.set_xlim(0, PANEL_W)
key_ax.set_ylim(0, 1)
key_ax.text(LEFT_IN - 0.30, 0.5, "Array length", ha="left", va="center")
x = LEFT_IN + 0.52
for L in (1e6, 3e6, 6e6):
    key_ax.scatter([x], [0.5], s=area(L), c=GREY, lw=0)
    key_ax.text(x + 0.07, 0.5, f"{L/1e6:g} Mb", ha="left", va="center")
    x += 0.42

fig.canvas.draw()
right_edge = max(t.get_window_extent(fig.canvas.get_renderer()).x1 for t in fig.texts + ax.texts + key_ax.texts) / fig.dpi
assert right_edge <= PANEL_W, f"labels run {right_edge - PANEL_W:.2f} in past the right edge"

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
