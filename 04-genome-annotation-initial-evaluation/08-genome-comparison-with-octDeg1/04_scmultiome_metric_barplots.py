#!/usr/bin/env python3
"""
Panel D of Figure 3 — single-cell Multiome metrics for the same four samples
mapped to OctDeg1.0 vs OctDeg2.0 (replaces the former separate panels F and G).
The OctDeg2.0 values come from the matched-annotation reruns (Liftoff-only reference;
07_cellranger_arc_count_octDeg2.sh, written into data/ by 08_update_fig3d_octdeg2_metrics.py):
  top:    feature linkages detected
  bottom: median genes per cell
Grouped bars (grey = OctDeg1.0, blue = OctDeg2.0) with one shared sample axis
and one assembly key. Drawn at printed size for Genome Research (8 pt text),
in the same format as panel C.

The samples of interest are the scMultiome samples in TARGET_SAMPLES; the full
feature_linkage table holds many more samples and is filtered to these.

Inputs (bundled in data/):
  feature_linkage.txt   samples / genome / "Feature linkages detected"
  gene_per_cell.txt     samples / genome / "Median genes per cell"

Run with the python-visualizations env:
    python 04_scmultiome_metric_barplots.py [--width-mm 84] [--height-mm 60]
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


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


FONT_SIZE = 8            # pt at final print size; Genome Research: 8-10 pt
TARGET_SAMPLES = ["181_PFC", "181_dHIP", "6997_dHIP", "7000_PFC"]
ASSEMBLIES = [("octDeg1", "OctDeg1.0", _cfg("FIG3_COLOR_OCTDEG1")),
              ("octDeg2", "OctDeg2.0", _cfg("FIG3_COLOR_OCTDEG2"))]
SAMPLE_LABEL = lambda s: s.replace("_", " ")        # display name: no underscores

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=84)
ap.add_argument("--height-mm", type=float, default=60)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4

DATA = Path(__file__).resolve().parent / "data"
out_dir = _cfg("FIG3_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, "figure3D_multiome_metrics")


def load(fname, value_col):
    d = pd.read_csv(DATA / fname, sep="\t").drop_duplicates()
    d = d[d["samples"].isin(TARGET_SAMPLES)]
    t = d.pivot(index="samples", columns="genome", values=value_col).loc[TARGET_SAMPLES]
    assert t.notna().all().all(), f"{fname}: missing values for {TARGET_SAMPLES}"
    return t


METRICS = [  # (table, y-axis label, scale divisor, y ticks)
    (load("feature_linkage.txt", "Feature linkages detected"), "Feature linkages\n(millions)", 1e6,
     [0, 0.5, 1, 1.5]),
    (load("gene_per_cell.txt", "Median genes per cell"), "Median genes\nper cell (thousands)", 1e3,
     [0, 1, 2, 3]),
]
for t, *_ in METRICS:
    print(t)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN, GAP_IN = 0.55, 0.05, 0.22, 0.20, 0.10
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
h = (PANEL_H - TOP_IN - BOTTOM_IN - GAP_IN) / len(METRICS)
x = np.arange(len(TARGET_SAMPLES))
width = 0.36
for k, (tbl, ylabel, div, yt) in enumerate(METRICS):
    y0 = PANEL_H - TOP_IN - (k + 1) * h - k * GAP_IN
    ax = fig.add_axes((LEFT_IN / PANEL_W, y0 / PANEL_H, 1 - (LEFT_IN + RIGHT_IN) / PANEL_W, h / PANEL_H))
    for i, (key, _, color) in enumerate(ASSEMBLIES):
        ax.bar(x + (i - 0.5) * width, tbl[key] / div, width, color=color, lw=0, zorder=3)
    top = max(yt[-1], ax.dataLim.y1) * 1.06          # never clip a bar
    ax.set_ylim(0, top)
    ax.set_yticks(yt, [f"{v:g}" for v in yt])
    ax.set_ylabel(ylabel, labelpad=3, linespacing=1.0)
    ax.yaxis.set_label_coords(-0.30 / (PANEL_W - LEFT_IN - RIGHT_IN), 0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([SAMPLE_LABEL(s) for s in TARGET_SAMPLES] if k == len(METRICS) - 1 else [])
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_xlim(-0.6, len(TARGET_SAMPLES) - 0.4)
    ax.yaxis.grid(True, color="#E5E5E5", lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

fig.legend(handles=[Patch(color=c, label=l) for _, l, c in ASSEMBLIES], ncol=2,
           loc="upper right", bbox_to_anchor=(1 - RIGHT_IN / PANEL_W, 1), frameon=False,
           handlelength=1.2, handleheight=0.9, columnspacing=1.2, borderaxespad=0.1)

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
