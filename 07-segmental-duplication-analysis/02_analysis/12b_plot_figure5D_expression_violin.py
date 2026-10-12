#!/usr/bin/env python3
"""
Figure 5D: expression of the eight GO:0071346 genes and their paralog copies across
29 RNA-seq samples from six degu tissues.

Same data and TPM calculation as 12_expression_violin_for_genes_of_interest.ipynb:
featureCounts gene counts (merged_paralog_singlecopy_counts.tsv) -> TPM within the counted
gene set, gene membership from merged_paralog_singlecopy_families.tsv (the same table used
for Figure 5C). One small plot per gene family, one violin per gene/paralog copy, points =
samples coloured by tissue.

Revised for Genome Research: 8 pt, italic gene symbols (titles and tick labels, not bold),
Okabe-Ito tissue colours, light-grey violins, one shared y-axis title and tissue key.

Run with the python-visualizations env:
    python 12b_plot_figure5D_expression_violin.py [--width-mm 175] [--height-mm 72]
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


FONT_SIZE = 8
ROWS = [["Gbp1", "Gbp5", "Gbp6"], ["Ifng", "Casp1", "Ccl3", "Lgals9", "Syncrip"]]
SAMPLE_TO_TISSUE = {
    "SRR17216293": "skin", "SRR17216294": "brain", "SRR17216295": "heart",
    "SRR17216296": "kidney", "SRR17216297": "liver", "SRR17216298": "lung",
    "SRR17216299": "skin", "SRR17216300": "brain", "SRR17216301": "heart",
    "SRR17216302": "kidney", "SRR17216303": "liver", "SRR17216304": "lung",
    "SRR17216305": "skin", "SRR17216306": "brain", "SRR17216307": "heart",
    "SRR17216308": "kidney", "SRR17216309": "liver", "SRR17216310": "lung",
    "SRR17216311": "brain", "SRR17216312": "heart", "SRR17216313": "kidney",
    "SRR17216314": "kidney", "SRR17216315": "lung", "SRR17216316": "skin",
    "SRR17216317": "heart", "SRR17216318": "kidney", "SRR17216319": "liver",
    "SRR17216320": "brain", "SRR17216321": "lung",
}
TISSUE_COLORS = {  # Okabe-Ito
    "brain": "#0072B2", "heart": "#D55E00", "kidney": "#009E73",
    "liver": "#E69F00", "lung": "#CC79A7", "skin": "#56B4E9",
}
VIOLIN_FC = "#E3E3E3"

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=175)
ap.add_argument("--height-mm", type=float, default=72)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4
here = os.path.dirname(os.path.abspath(__file__))
out_dir = _cfg("FIG5_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, "figure5D_expression_violin")

# ---------------------------------------------------------------- data (notebook TPM)
counts = pd.read_csv(os.path.join(here, "merged_paralog_singlecopy_counts.tsv"), sep="\t", comment="#", index_col=0)
samples = [c for c in counts.columns if c.startswith("SRR")]
fams = pd.read_csv(os.path.join(here, "merged_paralog_singlecopy_families.tsv"), sep="\t")
meta = fams[["gene_id", "gene_name", "length", "family"]].drop_duplicates().set_index("gene_id")
df = counts[samples].reset_index().melt(id_vars="Geneid", var_name="Sample", value_name="count")
df = df.rename(columns={"Geneid": "gene_id"})
df["length"] = df["gene_id"].map(meta["length"])
df["RPK"] = df["count"] / (df["length"] / 1000.0)
df["TPM"] = df["RPK"] / df.groupby("Sample")["RPK"].transform("sum") * 1e6
df["Gene"] = df["gene_id"].map(meta["gene_name"])
df["family"] = df["gene_id"].map(meta["family"])
df = df.dropna(subset=["Gene"])
df["Sample"] = df["Sample"].str.replace("_Aligned.sortedByCoord.out.bam", "", regex=False)
df["Tissue"] = df["Sample"].map(SAMPLE_TO_TISSUE)
assert df["Tissue"].notna().all(), "sample without tissue"
print(f"samples: {df['Sample'].nunique()}, tissues: {df['Tissue'].value_counts().div(df['Gene'].nunique()).astype(int).to_dict()}")


def members(fam):
    g = sorted(df.loc[df["family"] == fam, "Gene"].unique(), key=lambda x: (x != fam, x))
    assert g, f"no expression rows for {fam}"
    return g


# ---------------------------------------------------------------- plot
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE, "axes.titlesize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 0.40, 0.62, 0.14, 0.02
ROW_GAP_IN, COL_GAP_IN, XLAB_IN = 0.10, 0.28, 0.50      # XLAB_IN = room for rotated tick labels per row
rng = np.random.default_rng(1)
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
row_h = (PANEL_H - TOP_IN - BOTTOM_IN - ROW_GAP_IN - 2 * XLAB_IN - 2 * 0.13) / 2   # 0.13 = title
avail_w = PANEL_W - LEFT_IN - RIGHT_IN
max_slots = max(sum(len(members(f)) for f in r) for r in ROWS)
slot_w = (avail_w - (max(len(r) for r in ROWS) - 1) * COL_GAP_IN) / max_slots

axes = []
for ri, row in enumerate(ROWS):
    y0 = PANEL_H - TOP_IN - 0.13 - (ri + 1) * row_h - ri * (ROW_GAP_IN + XLAB_IN + 0.13)
    x = LEFT_IN
    for fam in row:
        genes = members(fam)
        w = len(genes) * slot_w
        ax = fig.add_axes((x / PANEL_W, y0 / PANEL_H, w / PANEL_W, row_h / PANEL_H))
        data = [df.loc[df["Gene"] == g, "TPM"].to_numpy() for g in genes]
        pos = np.arange(len(genes))
        varied = [i for i, d in enumerate(data) if d.max() > d.min()]    # KDE needs >1 distinct value
        if varied:
            parts = ax.violinplot([data[i] for i in varied], positions=pos[varied], widths=0.8,
                                  showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(VIOLIN_FC); b.set_edgecolor("#9A9A9A"); b.set_linewidth(0.4); b.set_alpha(1)
        for i, g in enumerate(genes):
            sub = df[df["Gene"] == g]
            ax.scatter(i + rng.uniform(-0.2, 0.2, len(sub)), sub["TPM"], s=4,
                       c=sub["Tissue"].map(TISSUE_COLORS), lw=0, zorder=3)
        ax.set_xticks(pos)
        ax.set_xticklabels(genes, rotation=45, ha="right", rotation_mode="anchor", fontstyle="italic")
        ax.set_xlim(-0.6, len(genes) - 0.4)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3, min_n_ticks=3))
        ax.set_title(fam, fontstyle="italic", pad=2)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        axes.append(ax)
        x += w + COL_GAP_IN

ytitle = fig.text(0.02 / PANEL_W, 0.5, "Expression (TPM)", rotation=90, ha="left", va="center")
handles = [Line2D([], [], ls="none", marker="o", ms=3.5, mfc=c, mec="none", label=t)
           for t, c in TISSUE_COLORS.items()]
fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.5), frameon=False,
           title="Tissue", alignment="left", handletextpad=0.2, labelspacing=0.35, borderaxespad=0.05)

fig.canvas.draw()
rend = fig.canvas.get_renderer()
ytitle_right = ytitle.get_window_extent(rend).x1
for ax in axes:
    assert all(t.get_window_extent(rend).x0 > ytitle_right for t in ax.get_yticklabels() if t.get_text()), \
        "y tick labels collide with the shared y-axis title"
    ylab_left = min(t.get_window_extent(rend).x0 for t in ax.get_yticklabels() if t.get_text())
    assert ylab_left >= 0, "y tick labels off the left edge"
    bot = min(t.get_window_extent(rend).y0 for t in ax.get_xticklabels())
    assert bot >= 0, "x tick labels off the bottom"
for a, b in zip(axes, axes[1:]):   # neighbouring plots in a row must not touch
    if abs(a.get_position().y0 - b.get_position().y0) < 1e-6:
        assert a.get_window_extent(rend).x1 < min(t.get_window_extent(rend).x0 for t in b.get_yticklabels() if t.get_text()), \
            "y tick labels collide with the plot on the left"

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
