#!/usr/bin/env python
"""Figure 2B, step 2: plot gene and repeat coverage along each chromosome.

Reads the cache from 02_prepare_figure2B_density.py. Drawn at printed size for
Genome Research: 8 pt text, chromosome names as one row of numbers (leader lines
where small chromosomes crowd), a single 50 Mb scale bar instead of per-chromosome
Mb ticks, no title. All three tracks share a 0-100% axis.

Run with the python-visualizations env:
    python 02_plot_figure2B_density.py [--width-mm 115]
"""

import argparse
import os
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec


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


FONT_SIZE = 8           # pt at final print size; Genome Research: 8-10 pt
PANEL_H = 2.3           # in (~58 mm)
CHROMS = [f"chr{i}" for i in range(1, 29)] + ["chrX", "chrY"]
TRACKS = [  # (track in the TSV, row label, color)
    ("Liftoff", "Liftoff", "#009E73"),   # Okabe-Ito bluish green
    ("BRAKER3", "BRAKER3", "#0072B2"),   # Okabe-Ito blue
    ("Repeats", "Repeats", "#666666"),
]
SCALE_MB = 100
LABEL_GAP_IN = 0.035    # in; min gap between neighbouring chromosome labels (~0.9 mm)
LEADER_MIN_IN = 0.015   # in; draw a leader line when a label moved more than this

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=115,
                help="panel width in mm (175 full width minus panel A and the 7 mm gap)")
ap.add_argument("--labels", choices=("vertical", "rows"), default="rows",
                help="chromosome numbers: vertical = one row rotated 90 deg; rows = two staggered rows")
opts = ap.parse_args()
PANEL_W = opts.width_mm / 25.4

out_dir = _cfg("FIG2B_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, f"figure2B_gene_repeat_density_{opts.width_mm:.0f}mm_{opts.labels}")

df = pd.read_csv(_cfg("FIG2B_DENSITY_TSV"), sep="\t")
lengths = df.groupby("chrom")["end"].max().reindex(CHROMS)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
    "axes.linewidth": 0.5, "ytick.major.width": 0.5, "ytick.major.size": 2,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

# Margins in inches, converted to figure fractions
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 0.62, 0.03, 0.05, 0.50
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
gs = gridspec.GridSpec(len(TRACKS), len(CHROMS), figure=fig, width_ratios=lengths.values,
                       left=LEFT_IN / PANEL_W, right=1 - RIGHT_IN / PANEL_W,
                       top=1 - TOP_IN / PANEL_H, bottom=BOTTOM_IN / PANEL_H,
                       wspace=0.25, hspace=0.18)

axes = [[None] * len(CHROMS) for _ in TRACKS]
for r, (track, label, color) in enumerate(TRACKS):
    sub = df[df["track"] == track]
    for c, chrom in enumerate(CHROMS):
        ax = fig.add_subplot(gs[r, c])
        axes[r][c] = ax
        d = sub[sub["chrom"] == chrom]
        ax.plot((d["start"] + d["end"]) / 2, d["pct_smoothed"], color=color, lw=0.6)
        ax.set_xlim(0, lengths[chrom])
        ax.set_ylim(0, 100)
        ax.set_xticks([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#999999")
        if c == 0:
            ax.set_yticks([0, 50, 100])
            ax.set_ylabel(label, labelpad=2)
        else:
            ax.set_yticks([])
            ax.spines["left"].set_visible(False)
        ax.yaxis.grid(True, color="#E5E5E5", lw=0.4)
        ax.set_axisbelow(True)

# Shared y-axis title, left of the row labels
fig.text(0.01 / PANEL_W, (BOTTOM_IN + (PANEL_H - TOP_IN - BOTTOM_IN) / 2) / PANEL_H,
         "Coverage per 1 Mb (%)", rotation=90, ha="left", va="center")

# Chromosome labels, directly under each chromosome (no leader lines):
#   vertical - one row, rotated 90 deg (a 2-digit 8 pt number is ~3 mm wide flat but ~2 mm on its side)
#   rows     - flat, alternating between two rows
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
bottom_axes = axes[-1]
y_axes_bottom = bottom_axes[0].get_position().y0
texts = []
for c, chrom in enumerate(CHROMS):
    pos = bottom_axes[c].get_position()
    x = (pos.x0 + pos.x1) / 2
    name = chrom.replace("chr", "")
    if opts.labels == "vertical":
        t = fig.text(x, y_axes_bottom - 0.05 / PANEL_H, name, ha="center", va="top", rotation=90)
    else:
        row = c % 2
        t = fig.text(x, y_axes_bottom - (0.05 + 0.14 * row) / PANEL_H, name, ha="center", va="top")
    texts.append(t)
fig.canvas.draw()
boxes = [t.get_window_extent(renderer) for t in texts]
for b1, b2 in zip(boxes, boxes[1:] if opts.labels == "vertical" else boxes[2:]):
    assert b2.x0 >= b1.x1, "chromosome labels overlap; widen the panel"
label_bottom = min(t.get_window_extent(renderer).y0 for t in texts) / fig.dpi / PANEL_H

# "Chromosome" title and a 50 Mb scale bar (same bp-to-length scale as every chromosome)
x_title = (bottom_axes[0].get_position().x0 + bottom_axes[-1].get_position().x1) / 2
fig.text(x_title, label_bottom - 0.05 / PANEL_H, "Chromosome", ha="center", va="top")
ax0 = bottom_axes[0]
frac_per_bp = (ax0.get_position().x1 - ax0.get_position().x0) / lengths[CHROMS[0]]
bar_x1 = bottom_axes[-1].get_position().x1
bar_x0 = bar_x1 - SCALE_MB * 1e6 * frac_per_bp
bar_y = label_bottom - 0.105 / PANEL_H   # same line as the "Chromosome" title
fig.add_artist(plt.Line2D([bar_x0, bar_x1], [bar_y, bar_y], transform=fig.transFigure, color="black", lw=1))
fig.text(bar_x0 - 0.03 / PANEL_W, bar_y, f"{SCALE_MB} Mb", ha="right", va="center")

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(out_stem + ".{pdf,svg,png}")
