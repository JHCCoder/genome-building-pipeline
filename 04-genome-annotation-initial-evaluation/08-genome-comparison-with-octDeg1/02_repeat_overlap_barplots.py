#!/usr/bin/env python3
"""
Panels B and C of Figure 3 — how the OctDeg2.0 -> OctDeg1.0 alignment behaves
with respect to the repetitive / non-repetitive regions of OctDeg2.0.

  Panel B: of the ALIGNED basepairs, the % that fall within repetitive vs
           non-repetitive regions of OctDeg2.0.
  Panel C: of the repetitive vs non-repetitive REGIONS of OctDeg2.0, the %
           that are aligned to OctDeg1.0.

Computation (equivalent to figure/barplots-feature-overlaps-with-repeat/
alignment_in_repeat_region.ipynb, reimplemented in pure pandas/numpy so no
bedtools/pyranges is required):
  1. merge the aligned query-side intervals,
  2. merge the RepeatMasker intervals,
  3. exact basepair overlap via a two-pointer sweep over the merged intervals,
  4. express the overlap against two denominators (total aligned bp for
     Panel B; total repetitive / non-repetitive genome bp for Panel C).

Inputs (see PATHS below):
  * aligned query-side BED (chrom, start, end) — OctDeg2.0 side
  * RepeatMasker .out.gff (chrom, start, end) for the final OctDeg2.0 assembly
  * chromosome sizes (.fai) for the final OctDeg2.0 assembly

Outputs (written next to this script):
  figure3B_aligned_in_repeat.png
  figure3C_repeat_aligned.png

Run with an env that has pandas/numpy/matplotlib (e.g. genome-assembly):

    conda activate genome-assembly
    python 02_repeat_overlap_barplots.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
# Project root — edit for your environment.
PROJ_ROOT = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"

# This script's directory (all outputs are written here).
HERE = Path(__file__).resolve().parent

ALIGNED_BED = (
    f"{PROJ_ROOT}/figure/align-octDeg1-hifiasm-minimap2/hifi041425_alignment.bed"
)
REPEAT_GFF = (
    f"{PROJ_ROOT}/figure/circos-plot/feature-overview/"
    "assembly_final.sorted.headerRenamed.fasta.out.chr.gff"
)
CHROM_FAI = (
    f"{PROJ_ROOT}/data/denovo_OctDegus_genome/041425-assembly/"
    "hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/"
    "assembly_final.sorted.headerRenamed.chrAssigned.fasta.fai"
)

OUT_B = str(HERE / "figure3B_aligned_in_repeat.png")
OUT_C = str(HERE / "figure3C_repeat_aligned.png")

# ============================================================
# Interval utilities (0-based, half-open [start, end))
# ============================================================
def merge_intervals(df):
    """Merge overlapping/touching intervals within each chromosome."""
    df = df[["chrom", "start", "end"]].sort_values(["chrom", "start", "end"])
    rows = []
    for chrom, g in df.groupby("chrom", sort=False):
        s = g["start"].values
        e = g["end"].values
        cur_s, cur_e = s[0], e[0]
        for i in range(1, len(s)):
            if s[i] <= cur_e:                 # overlap or touching
                cur_e = max(cur_e, e[i])
            else:
                rows.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s[i], e[i]
        rows.append((chrom, cur_s, cur_e))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def overlap_bp(a, b):
    """Total basepair overlap between two merged (non-overlapping, sorted)
    interval sets, matched by chromosome, via a two-pointer sweep."""
    b_by_chrom = {c: g for c, g in b.groupby("chrom", sort=False)}
    total = 0
    for chrom, ga in a.groupby("chrom", sort=False):
        if chrom not in b_by_chrom:
            continue
        gb = b_by_chrom[chrom]
        as_ = ga["start"].values
        ae = ga["end"].values
        bs_ = gb["start"].values
        be = gb["end"].values
        i = j = 0
        while i < len(as_) and j < len(bs_):
            lo = max(as_[i], bs_[j])
            hi = min(ae[i], be[j])
            if lo < hi:
                total += hi - lo
            if ae[i] < be[j]:
                i += 1
            else:
                j += 1
    return int(total)


# ============================================================
# Load + merge intervals
# ============================================================
aligned = pd.read_csv(ALIGNED_BED, sep="\t", names=["chrom", "start", "end"])
merged_aligned = merge_intervals(aligned)

repeat = pd.read_csv(
    REPEAT_GFF, sep="\t", header=None,
    names=["chrom", "source", "type", "start", "end",
           "score", "strand", "phase", "attributes"],
)[["chrom", "start", "end"]]
merged_repeat = merge_intervals(repeat)

# ============================================================
# Basepair overlap + denominators
# ============================================================
overlap_bp_ = overlap_bp(merged_aligned, merged_repeat)
total_aligned_bp = int((merged_aligned["end"] - merged_aligned["start"]).sum())
total_repeat_bp = int((merged_repeat["end"] - merged_repeat["start"]).sum())
non_overlap_bp = total_aligned_bp - overlap_bp_

# Total (final-assembly) genome size restricted to chromosomes present in the
# RepeatMasker output, so the non-repetitive denominator is consistent.
fai = pd.read_csv(CHROM_FAI, sep="\t", header=None,
                  names=["chrom", "length", "x", "y", "z"])
chroms_with_repeat = merged_repeat["chrom"].unique()
total_genome_bp = int(fai[fai["chrom"].isin(chroms_with_repeat)]["length"].sum())
total_nonrepeat_bp = total_genome_bp - total_repeat_bp

print(f"total aligned bp         : {total_aligned_bp:,}")
print(f"total repetitive bp      : {total_repeat_bp:,}")
print(f"total non-repetitive bp  : {total_nonrepeat_bp:,}")
print(f"overlap bp               : {overlap_bp_:,}")
print(f"aligned in repetitive    : {overlap_bp_/total_aligned_bp*100:.2f}%")
print(f"aligned in non-repetitive: {non_overlap_bp/total_aligned_bp*100:.2f}%")
print(f"repetitive aligned       : {overlap_bp_/total_repeat_bp*100:.2f}%")
print(f"non-repetitive aligned   : {non_overlap_bp/total_nonrepeat_bp*100:.2f}%")

# ============================================================
# Barplots
# ============================================================
def two_bar_panel(rep_value, nonrep_value, ylabel, out_path):
    fig, ax = plt.subplots(figsize=(5.5, 6.5))
    ax.bar(0, nonrep_value, 0.3, color="lightblue", alpha=0.85)
    ax.bar(0.4, rep_value, 0.3, color="lightcoral", alpha=0.85)
    for x, v in [(0, nonrep_value), (0.4, rep_value)]:
        ax.text(x, v + 1, f"{v:.1f}%", ha="center", va="bottom", fontsize=14)
    ax.text(0, -1, "Non-repetitive\nregion", ha="center", va="top", fontsize=12,
            color="darkblue", weight="bold")
    ax.text(0.4, -1, "Repetitive\nregion", ha="center", va="top", fontsize=12,
            color="darkred", weight="bold")
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_xticks([])
    ax.set_ylim(0, max(rep_value, nonrep_value) * 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.15, 0.55)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("saved:", out_path)


# Panel B — % of aligned regions within rep / non-rep
two_bar_panel(
    overlap_bp_ / total_aligned_bp * 100,
    non_overlap_bp / total_aligned_bp * 100,
    "% of OctDeg2.0 aligned regions\nin repetitive vs non-repetitive regions",
    OUT_B,
)

# Panel C — % of rep / non-rep regions aligned to OctDeg1.0
two_bar_panel(
    overlap_bp_ / total_repeat_bp * 100,
    non_overlap_bp / total_nonrepeat_bp * 100,
    "% of repetitive vs non-repetitive\nOctDeg2.0 regions aligned to OctDeg1.0",
    OUT_C,
)
