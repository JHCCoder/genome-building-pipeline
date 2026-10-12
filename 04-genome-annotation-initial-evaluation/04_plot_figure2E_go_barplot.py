#!/usr/bin/env python
"""Figure 2E: overlap fraction of the leading GO BP terms of selected clusters from 2D.

Reads the cache from 04_prepare_figure2DE_go_clusters.R. The terms in TERMS (from the
steroid metabolic process cluster and the clusters ranked from antimicrobial humoral
response down to MHC class II; chosen with the author on 2026-10-11) are drawn as bars of
overlap fraction (novel genes in the term / term size), filled with their cluster's colour
from 2D and grouped by cluster. Drawn at printed size for
Genome Research (8 pt text, no title); labels wrap to fit the label column.

Run with the python-visualizations env:
    python 04_plot_figure2E_go_barplot.py [--width-mm 38]
"""

import argparse
import os
import subprocess
import sys
import textwrap
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.transforms import Bbox


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
LINE_IN = FONT_SIZE * 1.15 / 72     # one wrapped label line
ROW_PAD_IN = 0.015      # extra space per row
GROUP_GAP_IN = 0.05     # extra space between clusters
LABEL_FRAC = 0.62       # share of the panel width used by the label column; bars use the rest
BAR_H_IN = 0.085        # bar thickness
TERMS = [  # shown in this order, grouped by their 2D cluster (chosen with the author, 2026-10-11)
    "steroid metabolic process",               # 1 steroid metabolic process
    "sterol metabolic process",
    "steroid hormone biosynthetic process",
    "antimicrobial humoral response",          # 6 antimicrobial humoral response
    "defense response to virus",
    "epoxygenase P450 pathway",                # 7 epoxygenase P450 pathway
    "arachidonate metabolic process",
    "inflammatory response",                   # 8 inflammatory response
    "MHC class II protein complex assembly",   # 9 MHC class II protein complex assembly
]

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=38)
opts = ap.parse_args()
PANEL_W = opts.width_mm / 25.4
LABEL_W_IN = PANEL_W * LABEL_FRAC

out_dir = _cfg("FIG2DE_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, f"figure2E_go_barplot_{opts.width_mm:.0f}mm")

d = pd.read_csv(_cfg("FIG2DE_CLUSTERS_TSV"), sep="\t")
missing = set(TERMS) - set(d["term"])
assert not missing, f"not in the significant-term cache: {missing}"
d = d.set_index("term").loc[TERMS].reset_index()
assert d["top_rank"].notna().all(), "every selected term must belong to a cluster coloured in 2D"
print(d[["top_rank", "parentTerm", "term", "Overlap", "score"]].round(1).to_string(index=False))

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE, "xtick.labelsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "xtick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

# Wrap each label to the label column (measured with the real font)
probe = plt.figure(figsize=(PANEL_W, 1))
renderer = probe.canvas.get_renderer()


def text_w_in(s):
    t = probe.text(0, 0, s, fontsize=FONT_SIZE)
    w = t.get_window_extent(renderer).width / probe.dpi
    t.remove()
    return w


def wrap(s, max_in):
    for width in range(60, 8, -1):
        lines = textwrap.wrap(s, width)
        if all(text_w_in(l) <= max_in for l in lines):
            return "\n".join(lines)
    return s


d["label"] = [wrap(t, LABEL_W_IN - 0.06) for t in d["term"]]
plt.close(probe)

# Row positions, top to bottom, in inches from the top of the plot area
heights = [LINE_IN * (lab.count("\n") + 1) + ROW_PAD_IN for lab in d["label"]]
y, centres = 0.0, []
for i, h in enumerate(heights):
    if i and d.loc[i, "top_rank"] != d.loc[i - 1, "top_rank"]:
        y += GROUP_GAP_IN
    centres.append(y + h / 2)
    y += h
plot_h = y
AXIS_IN = 0.42                       # x tick labels + axis title
PANEL_H = plot_h + AXIS_IN + 0.04

fig = plt.figure(figsize=(PANEL_W, PANEL_H))
ax = fig.add_axes((LABEL_W_IN / PANEL_W, AXIS_IN / PANEL_H,
                   (PANEL_W - LABEL_W_IN - 0.10) / PANEL_W, plot_h / PANEL_H))
ax.set_ylim(plot_h, 0)
ax.set_xlim(0, 1)
ax.barh(centres, d["overlap_fraction"], height=BAR_H_IN, color=d["color"], lw=0, zorder=3)
ax.set_yticks(centres)
ax.set_yticklabels(d["label"], linespacing=1.0)
ax.tick_params(axis="y", length=0, pad=3)
ax.set_xticks([0, 0.5, 1])
ax.set_xticklabels(["0", "0.5", "1"])
ax.set_xlabel("Overlap\nfraction", labelpad=2, linespacing=1.0)
ax.xaxis.grid(True, color="#E5E5E5", lw=0.5, zorder=0)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)

fig.canvas.draw()
r = fig.canvas.get_renderer()
left_edge = min(t.get_window_extent(r).x0 for t in ax.get_yticklabels()) / fig.dpi
assert left_edge >= 0, f"labels run {-left_edge:.2f} in past the left edge"
print(f"panel: {opts.width_mm:.0f} x {PANEL_H * 25.4:.1f} mm ({PANEL_W:.2f} x {PANEL_H:.2f} in), "
      f"{len(d)} bars")

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(out_stem + ".{pdf,svg,png}")
