#!/usr/bin/env python3
"""
Panels F and G of Figure 3 — scMultiome (10x Multiome) alignment metrics for
the same four samples mapped to OctDeg1.0 vs OctDeg2.0.

  Panel F: feature linkages (Cell Ranger ARC "Feature linkages detected")
  Panel G: median genes per cell

Both are grouped barplots (one bar per assembly per sample). The four
samples-of-interest are the scMultiome samples (181_PFC, 181_dHIP, 6997_dHIP,
7000_PFC); the full feature_linkage table holds many more samples and is
filtered down to these four here.

Inputs (bundled in data/ for portability; canonical copies live at
figure/sample_cCRE_cisTrans_TAD_metrics/):

  feature_linkage.txt   samples / genome / "Feature linkages detected"
  gene_per_cell.txt     samples / genome / "Median genes per cell"

Outputs (written next to this script):
  figure3F_feature_linkage.png
  figure3G_genes_per_cell.png

Port of figure/sample_cCRE_cisTrans_TAD_metrics/
scatterplot_feature_linkage_scMultiome.ipynb (barplot cells).

Run with an env that has pandas/matplotlib:

    conda activate genome-assembly
    python 04_scmultiome_metric_barplots.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# RUN PARAMETERS
# ---------------------------------------------------------------------------
TARGET_SAMPLES = ["181_PFC", "181_dHIP", "6997_dHIP", "7000_PFC"]
ASSEMBLY_COLORS = {"octDeg1": "#33ffff", "octDeg2": "#ffcc00"}
DPI = 300

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

OUT_F = str(HERE / "figure3F_feature_linkage.png")
OUT_G = str(HERE / "figure3G_genes_per_cell.png")


def grouped_barplot(df, value_col, ylabel, out_path, fmt=None):
    """Grouped barplot of `value_col` per sample (one bar per assembly)."""
    df = df.sort_values(["samples", "genome"]).reset_index(drop=True)
    samples = df["samples"].unique()
    genomes = df["genome"].unique()

    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(len(samples))
    width = 0.35
    for i, g in enumerate(genomes):
        vals = df[df["genome"] == g][value_col]
        ax.bar(x + i * width, vals, width, label=g,
               color=ASSEMBLY_COLORS[g], alpha=0.7)
    ax.set_xlabel("Samples", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(samples, fontsize=14)
    ax.tick_params(axis="y", labelsize=14)
    if fmt == "millions":
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, p: f"{v * 1e-6:.0f}M"))
    elif fmt == "thousands":
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, p: f"{v * 1e-3:.0f}K"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=ASSEMBLY_COLORS[g],
                             alpha=0.7, label=g) for g in genomes]
    ax.legend(handles=handles, title="Assembly", fontsize=13, title_fontsize=14,
              loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out_path)


# ============================================================
# Panel F — feature linkages (filter to the 4 samples of interest)
# ============================================================
fl = pd.read_csv(DATA / "feature_linkage.txt", sep="\t", header=0)
fl = fl[fl["samples"].isin(TARGET_SAMPLES)]
print("feature-linkage samples:", list(fl["samples"].unique()))
grouped_barplot(fl, "Feature linkages detected", "Feature linkages",
                OUT_F, fmt="millions")

# ============================================================
# Panel G — median genes per cell
# ============================================================
gpc = pd.read_csv(DATA / "gene_per_cell.txt", sep="\t", header=0)
print("genes-per-cell samples:", list(gpc["samples"].unique()))
grouped_barplot(gpc, "Median genes per cell", "Median genes per cell",
                OUT_G, fmt="thousands")
