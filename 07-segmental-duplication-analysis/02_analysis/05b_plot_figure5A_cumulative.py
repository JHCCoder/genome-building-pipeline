#!/usr/bin/env python3
"""
Figure 5A: cumulative percentage of SegDup length and of SegDup-overlapping features,
with merged SegDup regions ordered from longest to shortest.

Same quantities as the final plotting cell of 05_cumulative_feature_overlap.ipynb, read
from that notebook's cached per-region counts ($FIG5_SEGDUP_COUNTS_TSV; 28,961 merged
regions on assembled chromosomes):
  * SegDup length (black) — the 80% point defines the SegDup hotspots;
  * unique genes covered (each gene counted once, the first time any region covers it);
  * repeats (all, LINE, SINE, DNA transposon, unclassified).
Hotspots = the regions whose cumulative length is <= 80% (1,425 regions), the same set used
for the GO enrichment (B) and the circos (C); the marker and labels are placed there.

Drawn at printed size for Genome Research (8 pt, no title, Okabe-Ito colours).

Run with the python-visualizations env:
    python 05b_plot_figure5A_cumulative.py [--width-mm 82] [--height-mm 36]
"""

import argparse
import os
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from intervaltree import IntervalTree


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


FONT_SIZE = 8
CURVES = [  # (column, label, colour, line width)
    ("Length", "SegDup length", "#000000", 1.4),
    ("gene", "Genes", "#0072B2", 1.2),
    ("repeat_count", "All repeats", "#D55E00", 1.0),
    ("LINE_repeat_count", "LINE", "#E69F00", 0.9),
    ("SINE_repeat_count", "SINE", "#CC79A7", 0.9),
    ("DNA_transposon_repeat_count", "DNA transposon", "#009E73", 0.9),
    ("unknown_repeat_count", "Unclassified", "#56B4E9", 0.9),
]
THRESHOLD = 80

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=82)
ap.add_argument("--height-mm", type=float, default=36)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4
out_dir = _cfg("FIG5_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, "figure5A_cumulative_segdup")

seg = pd.read_csv(_cfg("FIG5_SEGDUP_COUNTS_TSV"), sep="\t")
seg = seg.sort_values("Length", ascending=False, kind="stable").reset_index(drop=True)
x = np.arange(1, len(seg) + 1)
pct = {}
for col, *_ in CURVES:
    if col == "gene":
        continue
    pct[col] = seg[col].cumsum().to_numpy() / seg[col].sum() * 100

# unique genes covered, walking regions from longest to shortest
genes = pd.read_csv(_cfg("SEGDUP_GENE_OVERLAP_TSV"), sep="\t")
trees = {c: IntervalTree.from_tuples(zip(g["Start"].astype(int), g["End"].astype(int), g["gene"]))
         for c, g in genes.groupby("Chromosome")}
covered, cum = set(), []
for c, s, e in zip(seg["chrom"], seg["start"], seg["end"]):
    t = trees.get(c)
    if t is not None:
        covered.update(iv.data for iv in t.overlap(int(s), int(e)))
    cum.append(len(covered))
cum = np.array(cum)
pct["gene"] = cum / cum[-1] * 100

# Hotspots = regions whose cumulative length is <= 80% (the set used for GO enrichment
# and the circos, segdup_hotspots_80percLength.tsv): 1,425 regions, 79.995% of length.
i80 = int(np.flatnonzero(pct["Length"] <= THRESHOLD)[-1])
n_regions, n_genes, gene_pct = i80 + 1, int(cum[i80]), float(pct["gene"][i80])
print(f"regions: {len(seg):,}; hotspots (<= {THRESHOLD}% of SegDup length): {n_regions:,} regions "
      f"({pct['Length'][i80]:.2f}%); unique genes there: {n_genes:,} of {cum[-1]:,} ({gene_pct:.1f}%)")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 0.40, 1.18, 0.05, 0.36   # key sits in the right margin
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
ax = fig.add_axes((LEFT_IN / PANEL_W, BOTTOM_IN / PANEL_H,
                   1 - (LEFT_IN + RIGHT_IN) / PANEL_W, 1 - (TOP_IN + BOTTOM_IN) / PANEL_H))
for col, label, color, lw in reversed(CURVES):           # SegDup length drawn last, on top
    ax.plot(x, pct[col], color=color, lw=lw, label=label, zorder=3 if col == "Length" else 2)
ax.axvline(n_regions, color="#7F7F7F", ls="--", lw=0.6, zorder=1)
ax.plot([n_regions], [pct["Length"][i80]], "o", color="black", ms=3.5, zorder=4)
ax.plot([n_regions], [gene_pct], "s", color=CURVES[1][2], ms=3.5, zorder=4)
ax.annotate(f"{n_regions:,} regions", (n_regions, pct["Length"][i80]), xytext=(4, -9),
            textcoords="offset points", ha="left", va="center")
ax.annotate(f"{n_genes:,} genes", (n_regions, gene_pct), xytext=(4, 0),
            textcoords="offset points", ha="left", va="center", color=CURVES[1][2])

ax.set_xlim(0, len(seg) * 1.01)
ax.set_ylim(0, 102)
ax.set_xticks([0, 10000, 20000], ["0", "10k", "20k"])
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("SegDup regions, longest to shortest", labelpad=2)
ax.set_ylabel("Cumulative (%)", labelpad=2)
ax.yaxis.grid(True, color="#E5E5E5", lw=0.5, zorder=0)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
          handlelength=1.2, labelspacing=0.35, borderaxespad=0)

fig.canvas.draw()
right = ax.get_legend().get_window_extent(fig.canvas.get_renderer()).x1 / fig.dpi
assert right <= PANEL_W, f"key runs {right - PANEL_W:.2f} in past the right edge"

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
