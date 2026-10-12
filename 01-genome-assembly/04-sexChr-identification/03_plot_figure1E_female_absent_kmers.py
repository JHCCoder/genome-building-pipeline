#!/usr/bin/env python
"""Figure 1E: percentage of each scaffold's k-mers absent from the female k-mer set.

DiscoverY output (03_kmer_overlap_analysis.sh), as plotted in cell 6 of
figure/chromosome-coverage-kmer-properties/chromosomal_rDiscover.ipynb, re-drawn at
printed size for Genome Research (8-10 pt text, no title, no in-panel note: the
definition belongs in the figure legend).

Scaffold N is the Nth longest of HapHiC scaffolds 1-30, the same numbering as
Figure 1D.

Run with the python-visualizations env:
    python 03_plot_figure1E_female_absent_kmers.py [--width column|full]
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
PANEL_H = 1.44          # in (~37 mm; 20% shorter than the first draft, to make room for panel F)
N_SCAFFOLDS = 30
BAR_COLOR = "#0072B2"   # Okabe-Ito blue

ap = argparse.ArgumentParser()
ap.add_argument("--width", choices=PANEL_WIDTHS, default="column",
                help="column = 84 mm (odd scaffolds labeled), full = 175 mm (all labeled)")
opts = ap.parse_args()
PANEL_W = PANEL_WIDTHS[opts.width]

out_dir = _cfg("FIG1E_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, f"figure1E_female_absent_kmers_{opts.width}")

kmer = pd.read_csv(_cfg("DISCOVERY_PROPORTIONS"), sep=" ",
                   names=["contig", "length", "proportion", "depth"])
kmer = kmer[kmer["contig"].str.match(r"PGA_scaffold_([1-9]|[12][0-9]|30)__")]
kmer = kmer.sort_values("length", ascending=False).reset_index(drop=True)
assert len(kmer) == N_SCAFFOLDS, f"expected {N_SCAFFOLDS} scaffolds, got {len(kmer)}"
kmer["scaffold"] = kmer.index + 1
kmer["percent"] = kmer["proportion"] * 100

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

fig, ax = plt.subplots(figsize=(PANEL_W, PANEL_H))
ax.bar(kmer["scaffold"], kmer["percent"], width=0.7, color=BAR_COLOR, lw=0, zorder=3)

ax.set_xlim(0.3, N_SCAFFOLDS + 0.7)
ax.set_ylim(0, 40)
ax.set_yticks(range(0, 41, 10))
ax.set_xticks(range(1, N_SCAFFOLDS + 1))
step = 2 if opts.width == "column" else 1
ax.set_xticklabels([str(i) if (i - 1) % step == 0 else "" for i in range(1, N_SCAFFOLDS + 1)])
ax.set_xlabel("Scaffold", labelpad=2)
ax.set_ylabel("k-mers absent from" + ("\n" if opts.width == "column" else " ") + "female (%)", labelpad=2)
ax.yaxis.grid(True, color="#E5E5E5", lw=0.5, zorder=0)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

fig.tight_layout(pad=0.15)

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(out_stem + ".{pdf,svg,png}")
