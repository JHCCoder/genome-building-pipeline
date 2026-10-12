#!/usr/bin/env python3
"""
Figure 5B: GO biological-process terms enriched among the SegDup-hotspot genes.

Reads the enrichment table written by 06_go_enrichment.ipynb (gseapy/Enrichr, GO BP 2025;
foreground = 469 non-LOC gene symbols in the 1,425 hotspot regions, paralog copies collapsed
to their parent symbol; background = 6,814 non-LOC symbols of SegDup-overlapping genes) and
plots every term with adjusted P < 0.05, ordered by combined score:
  x = Enrichr combined score, dot area = % of the term's genes in the hotspot set,
  colour = adjusted P (FDR), single-hue blue (darker = smaller FDR).
Replaces the gseapy dotplot (no title, no unexplained outer rings, no pale-yellow end of the
colour scale). Term names are shown in GO's lower-case style without IDs (IDs in the legend);
the term followed up in C and D is bold.

Run with the python-visualizations env:
    python 06b_plot_figure5B_go_dotplot.py [--width-mm 82] [--height-mm 60]
"""

import argparse
import os
import re
import subprocess
import sys
import textwrap
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize


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
FDR_MAX = 0.05
HIGHLIGHT = "GO:0071346"            # cellular response to type II interferon (followed up in C, D)
CMAP = LinearSegmentedColormap.from_list("fdr_blue", ["#08306B", "#2171B5", "#9ECAE1"])
SIZE_FOR = lambda pct: 4 + (pct / 100) ** 1.0 * 80     # dot area (pt^2) for a % of genes in term
WRAP = 40                                              # characters per label line

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=82)
ap.add_argument("--height-mm", type=float, default=60)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4
out_dir = _cfg("FIG5_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, "figure5B_go_dotplot")

here = os.path.dirname(os.path.abspath(__file__))
go = pd.read_csv(os.path.join(here, "go_enrichment_results_segdup_hotspot_seqRemoved_80percLength.csv"))
go = go[go["Adjusted P-value"] < FDR_MAX].copy()
go["go_id"] = go["Term"].str.extract(r"\((GO:\d+)\)", expand=False)
go["name"] = (go["Term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True).str.lower()
              .str.replace(r"\bii\b", "II", regex=True)
              .str.replace("natural killer", "NK", regex=False)          # standard abbreviations,
              .str.replace("interleukin-1 beta", "IL-1\u03b2", regex=False))  # full names in the legend
k_n = go["Overlap"].str.split("/", expand=True).astype(int)
go["pct"] = k_n[0] / k_n[1] * 100
go = go.sort_values(["Combined Score", "Adjusted P-value"], ascending=[True, False]).reset_index(drop=True)
print(go[["go_id", "name", "Overlap", "Adjusted P-value", "Combined Score"]].iloc[::-1].to_string(index=False))

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 0,
    "mathtext.fontset": "custom", "mathtext.rm": "Nimbus Sans", "mathtext.it": "Nimbus Sans:italic",
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 2.08, 0.10, 0.03, 0.74
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
ax = fig.add_axes((LEFT_IN / PANEL_W, BOTTOM_IN / PANEL_H,
                   1 - (LEFT_IN + RIGHT_IN) / PANEL_W, 1 - (TOP_IN + BOTTOM_IN) / PANEL_H))
norm = Normalize(0, FDR_MAX)
wrapped = [textwrap.fill(n, WRAP) for n in go["name"]]
row_h = np.array([1.0 + 0.45 * w.count("\n") for w in wrapped])   # two-line terms get a taller row
y = np.cumsum(row_h) - row_h / 2
ax.hlines(y, 0, go["Combined Score"], color="#E0E0E0", lw=0.5, zorder=1)
ax.scatter(go["Combined Score"], y, s=SIZE_FOR(go["pct"]), c=go["Adjusted P-value"], cmap=CMAP, norm=norm,
           edgecolors="none", zorder=3)
ax.set_yticks(y)
labels = ax.set_yticklabels(wrapped, linespacing=1.0)
for lab, gid in zip(labels, go["go_id"]):
    if gid == HIGHLIGHT:
        lab.set_fontweight("bold")
ax.set_ylim(-0.1, row_h.sum() + 0.1)
ax.set_xlim(0, go["Combined Score"].max() * 1.12)
ax.set_xlabel("Combined score", labelpad=2)
ax.xaxis.grid(True, color="#EEEEEE", lw=0.5, zorder=0)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

# Keys below the plot: dot size (% of term genes) and FDR colour bar
# Keys in one row below the x-axis title: dot size (left), FDR colour bar (right)
key = fig.add_axes((0.30 / PANEL_W, 0.02 / PANEL_H, 1.30 / PANEL_W, 0.30 / PANEL_H))
key.set_axis_off()
key.set_xlim(0, 1)
key.set_ylim(0, 1)
key.text(0.0, 0.98, "Term genes in set (%)", ha="left", va="top")
for xk, p in zip((0.06, 0.37, 0.68), (20, 50, 80)):
    key.scatter([xk], [0.24], s=SIZE_FOR(p), c="#7F7F7F", edgecolors="none")
    key.text(xk + 0.07, 0.24, f"{p}", ha="left", va="center")
BAR_X, BAR_W = 1.98, 0.95
cax = fig.add_axes((BAR_X / PANEL_W, 0.15 / PANEL_H, BAR_W / PANEL_W, 0.06 / PANEL_H))
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal",
                  ticks=[0, 0.025, 0.05])
cb.ax.set_xticklabels(["0", "0.025", "0.05"])
cb.outline.set_linewidth(0.5)
cb.ax.tick_params(width=0.5, length=2, pad=1)
cb.ax.get_xticklabels()[0].set_ha("left")    # end labels align with the bar ends
cb.ax.get_xticklabels()[-1].set_ha("right")
fig.text(BAR_X / PANEL_W, (0.15 + 0.06 + 0.02) / PANEL_H, "Adjusted $\\it{P}$ (FDR)", ha="left", va="bottom")

fig.canvas.draw()
rend = fig.canvas.get_renderer()
left = min(t.get_window_extent(rend).x0 for t in labels) / fig.dpi
assert left >= 0, f"term labels run {-left:.2f} in past the left edge"
right = max(t.get_window_extent(rend).x1 for t in fig.texts + key.texts) / fig.dpi
bottom = min(t.get_window_extent(rend).y0 for t in cb.ax.get_xticklabels()) / fig.dpi
assert bottom >= 0, "colour-bar tick labels run off the bottom"
xlab_bottom = ax.xaxis.label.get_window_extent(rend).y0
assert all(t.get_window_extent(rend).y1 < xlab_bottom for t in fig.texts + key.texts), "keys collide with the x-axis title"
boxes = sorted((t.get_window_extent(rend) for t in labels), key=lambda b: b.y0)
_gaps = [(b.y0 - a.y1) / fig.dpi * 25.4 for a, b in zip(boxes, boxes[1:])]
assert all(g >= -0.05 for g in _gaps), "term labels overlap"
assert right <= PANEL_W, f"keys run {right - PANEL_W:.2f} in past the right edge"

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
