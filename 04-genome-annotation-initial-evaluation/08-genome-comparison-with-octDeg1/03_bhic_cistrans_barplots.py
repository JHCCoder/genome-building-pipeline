#!/usr/bin/env python3
"""
Panels D and E of Figure 3 — bulk Hi-C (bHiC) alignment metrics for the same
four samples aligned to OctDeg1.0 vs OctDeg2.0.

  Panel D: cis contact count
  Panel E: cis/trans contact ratio

Both are grouped barplots (one bar per assembly per sample).

Input (bundled in data/ for portability; canonical copy lives at
figure/sample_cCRE_cisTrans_TAD_metrics/cis_trans_bhic_metrics.txt):

    Sample       Cis_Trans   Cis_Contact   Read_mapped   Pairs_mapped   Assembly
    degu_060302  1.237471031 2852961       12279254      5165880        octDeg1
    ...          ...         ...           ...           ...            ...

Outputs (written next to this script):
  figure3D_cis_contact.png
  figure3E_cistrans_ratio.png

Port of figure/sample_cCRE_cisTrans_TAD_metrics/barplot_cisTransRatio_bHiC.ipynb.

Run with an env that has pandas/matplotlib:

    conda activate genome-assembly
    python 03_bhic_cistrans_barplots.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# RUN PARAMETERS
# ---------------------------------------------------------------------------
ASSEMBLY_COLORS = {"octDeg1": "#33ffff", "octDeg2": "#ffcc00"}
DPI = 300

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "cis_trans_bhic_metrics.txt"

OUT_D = str(HERE / "figure3D_cis_contact.png")
OUT_E = str(HERE / "figure3E_cistrans_ratio.png")

# ============================================================
# Load data
# ============================================================
df = pd.read_csv(DATA, sep="\t", header=0)
samples = df["Sample"].unique()
assemblies = df["Assembly"].unique()
print("samples:", list(samples))
print("assemblies:", list(assemblies))


def grouped_barplot(metric, ylabel, out_path, millions=False):
    """Grouped barplot of `metric` per sample, one bar per assembly."""
    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(len(samples))
    width = 0.35
    for i, asm in enumerate(assemblies):
        vals = df[df["Assembly"] == asm][metric]
        ax.bar(x + i * width, vals, width, label=asm,
               color=ASSEMBLY_COLORS[asm], alpha=0.7)
    ax.set_xlabel("Samples", fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(samples, fontsize=14)
    ax.tick_params(axis="y", labelsize=14)
    if millions:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, p: f"{v * 1e-6:.0f}M"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=ASSEMBLY_COLORS[a],
                             alpha=0.7, label=a) for a in assemblies]
    ax.legend(handles=handles, title="Assembly", fontsize=13, title_fontsize=14,
              loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out_path)


# Panel D — cis contact
grouped_barplot("Cis_Contact", "Cis Contact", OUT_D, millions=True)

# Panel E — cis/trans ratio
grouped_barplot("Cis_Trans", "Cis/Trans Ratio", OUT_E)
