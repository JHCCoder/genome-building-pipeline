#!/usr/bin/env python
"""Figure 1D: mean read coverage per scaffold (1-30).

Long (HiFi) and short (Illumina) reads, each mapped to the full and to the
hard-masked assembly. Same data as cell 12 of 02_visualizing_averageChrCoverage.ipynb,
re-drawn at printed size for Genome Research (8-10 pt text, key inside the panel, no title).

Color = read type, fill = assembly (filled: full; open: hard-masked), so the
four series are told apart without relying on similar shades.

In all four tables, row N (the Nth longest scaffold) is scaffold N; the hard-masked
tables name it seqN, the full-assembly tables keep the HapHiC/Juicebox name.

Run with the python-visualizations env:
    python 02_plot_figure1D_coverage.py [--width column|full]
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
from matplotlib.lines import Line2D


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
PANEL_WIDTHS = {"column": 3.31, "full": 6.89}   # in; GR one column 84 mm, full width 175 mm
PANEL_H = 1.6           # in (~41 mm; 20% shorter than the first draft, to make room for panel F)
N_SCAFFOLDS = 30
MARKER_SIZE = 3.2       # pt

COL_LONG = "#0072B2"    # Okabe-Ito blue
COL_SHORT = "#D55E00"   # Okabe-Ito vermillion

STEM = "coverage_{read}_read_hifiasm_041425_scaffolded_juiceBox_sorted_{masked}primaryRead_contig.tsv"
SERIES = [  # (read, hard-masked?, label, color, marker)
    ("long", False, "Long reads", COL_LONG, "o"),
    ("long", True, "Long reads, masked", COL_LONG, "o"),
    ("short", False, "Short reads", COL_SHORT, "s"),
    ("short", True, "Short reads, masked", COL_SHORT, "s"),
]

ap = argparse.ArgumentParser()
ap.add_argument("--width", choices=PANEL_WIDTHS, default="column",
                help="column = 84 mm (odd scaffolds labeled), full = 175 mm (all labeled)")
opts = ap.parse_args()
PANEL_W = PANEL_WIDTHS[opts.width]

in_dir = _cfg("PRIMARY_COVERAGE_DIR")
out_dir = _cfg("FIG1D_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, f"figure1D_scaffold_coverage_{opts.width}")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
handles = []
for read, masked, label, color, marker in SERIES:
    path = os.path.join(in_dir, STEM.format(read=read, masked="hardMasked_" if masked else ""))
    cov = pd.read_csv(path, sep="\t", names=["contig", "start", "end", "coverage"]).head(N_SCAFFOLDS)
    x = range(1, len(cov) + 1)
    face = "white" if masked else color
    ax.plot(x, cov["coverage"], ls="none", marker=marker, ms=MARKER_SIZE, mfc=face, mec=color,
            mew=0.75, zorder=3)
    handles.append(Line2D([], [], ls="none", marker=marker, ms=MARKER_SIZE + 0.5, mfc=face,
                          mec=color, mew=0.75, label=label))

ax.set_xlim(0.3, N_SCAFFOLDS + 0.7)
ax.set_ylim(0, 190)   # headroom so the key clears scaffold 29 (133x, short reads, masked)
ax.set_xticks(range(1, N_SCAFFOLDS + 1))
step = 2 if opts.width == "column" else 1
ax.set_xticklabels([str(i) if (i - 1) % step == 0 else "" for i in range(1, N_SCAFFOLDS + 1)])
ax.set_yticks(range(0, 161, 40))
ax.set_xlabel("Scaffold", labelpad=2)
ax.set_ylabel("Mean coverage (×)", labelpad=2)
ax.yaxis.grid(True, color="#E5E5E5", lw=0.5, zorder=0)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=2, frameon=False,
          handletextpad=0.2, columnspacing=1.0, borderaxespad=0, labelspacing=0.3)

fig.tight_layout(pad=0.15)

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(out_stem + ".{pdf,svg,png}")
