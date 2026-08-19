#!/usr/bin/env python3
"""
Panel A of Figure 3 — whole-genome collinearity scatterplot between the new
OctDeg2.0 assembly (query; hifiasm-041425) and the OctDeg1.0 reference
(target), obtained with minimap2 (-x asm5).

Each point is one alignment. The x-axis concatenates the OctDeg2.0
chromosomes in natural order; the y-axis orders and orients the OctDeg1.0
contigs so that collinear alignments run along the diagonal. Points are
colored by the OctDeg2.0 chromosome they map to. Only high-confidence
alignments (mapq >= MAPQ_MIN, alignment length >= ALN_LEN_MIN) are shown, and
optionally only alignments matching each contig's dominant strand (syntenic).

Port of figure/align-octDeg1-hifiasm-minimap2/visualize_paf_clean_collinearity.ipynb.

Inputs (see PATHS below):
  * minimap2 -x asm5 PAF (query = OctDeg2.0, target = OctDeg1.0)

Output (written next to this script):
  figure3A_alignment_scatterplot.png

Run with an env that has pandas/numpy/natsort/matplotlib (e.g. genome-assembly):

    conda activate genome-assembly
    python 01_alignment_scatterplot.py
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ---------------------------------------------------------------------------
# RUN PARAMETERS
# ---------------------------------------------------------------------------
MAPQ_MIN = 60            # keep primary, high-confidence mappings only
ALN_LEN_MIN = 10_000     # bp; drop short / repetitive noise
SYNTENIC_ONLY = True     # keep only alignments matching each contig's dominant strand
DOT_SIZE = 3
DPI = 300

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
# Project root — edit for your environment.
PROJ_ROOT = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"

# This script's directory (all outputs are written here).
HERE = Path(__file__).resolve().parent

PAF = f"{PROJ_ROOT}/figure/align-octDeg1-hifiasm-minimap2/alignment_octDeg1_hifi041425.paf"

OUT_PNG = str(HERE / "figure3A_alignment_scatterplot.png")


def natural_key(s):
    """Natural sort key for chromosome names (chr1 < chr2 < ... < chr10 < chrX)."""
    return re.sub(r"(\d+)", lambda m: m.group(1).zfill(8), s)

# ============================================================
# Load + filter PAF
# ============================================================
cols = ["qname", "qlen", "qstart", "qend", "strand", "tname",
        "tlen", "tstart", "tend", "matches", "aln_len", "mapq"]
raw = pd.read_csv(PAF, sep="\t", header=None, usecols=range(12), names=cols)
raw = raw[~raw["qname"].str.contains("seq")].copy()          # drop tiny unplaced scaffolds
for c in ["qname", "tname"]:
    raw[c] = raw[c].astype(str)
for c in ["qstart", "qend", "tstart", "tend", "qlen", "tlen", "aln_len", "mapq", "matches"]:
    raw[c] = raw[c].astype(int)
raw["ident"] = raw["matches"] / raw["aln_len"]

df = raw[(raw["mapq"] >= MAPQ_MIN) & (raw["aln_len"] >= ALN_LEN_MIN)].copy()
print(f"raw alignments                     : {len(raw):,}")
print(f"after mapq>={MAPQ_MIN} & aln_len>={ALN_LEN_MIN:,}: {len(df):,} "
      f"({100*len(df)/len(raw):.1f}%)")
print(f"identity of kept alignments        : median={df['ident'].median():.3f}")

# ============================================================
# Contig assignment: primary chromosome + dominant strand
# ============================================================
grp = df.groupby(["tname", "qname", "strand"])["aln_len"].sum().rename("bp").reset_index()
dom_strand = (grp.sort_values(["tname", "bp"])
                 .drop_duplicates("tname", keep="last")
                 .set_index("tname")["strand"])
prim_chr = (df.groupby(["tname", "qname"])["aln_len"].sum()
               .sort_values().reset_index()
               .drop_duplicates("tname", keep="last")
               .set_index("tname")["qname"])
print(f"contigs with >=1 kept alignment: {len(prim_chr):,}")
print(f"reverse-dominant contigs        : {(dom_strand=='-').sum():,} "
      f"({100*(dom_strand=='-').mean():.0f}%)")

# ============================================================
# Axes: concatenate query (x); order + orient target contigs (y)
# ============================================================
chr_order = sorted(df["qname"].unique(), key=natural_key)
chr_len = df.groupby("qname")["qlen"].first().reindex(chr_order)
chr_x0 = chr_len.cumsum().shift(1, fill_value=0).astype(int)   # x offset of each chromosome
chr_rank = {c: i for i, c in enumerate(chr_order)}
GENOME_X = int(chr_len.sum())
print("chromosomes on x:", ", ".join(chr_order))

med = df.groupby(["tname", "qname"])["qstart"].median().rename("med_q").reset_index()
ord_tbl = med.merge(prim_chr.rename("pc"), left_on="tname", right_index=True)
ord_tbl = ord_tbl[ord_tbl["qname"] == ord_tbl["pc"]]             # position on primary chr
ord_tbl["rk"] = ord_tbl["pc"].map(chr_rank)
ord_tbl = ord_tbl.sort_values(["rk", "med_q"])
sorted_contigs = ord_tbl["tname"].tolist()

tlen = df.groupby("tname")["tlen"].first()
y_off, y_total = {}, 0
for t in sorted_contigs:
    y_off[t] = y_total
    y_total += int(tlen[t])
df["y_off"] = df["tname"].map(y_off).astype(int)
print(f"y-axis total: {y_total:,} bp across {len(sorted_contigs):,} contigs")

# Flip each contig so its dominant orientation reads 5'->3' with the hifiasm
# assembly. After this, collinear alignments all have positive slope.
dom = dom_strand[df["tname"]].values
s = df["strand"].values
ts, te = df["tstart"].values, df["tend"].values
tl = df["tlen"].values
fq = np.where(s == "+", ts, te)      # target-forward coordinate aligned to qstart
fqe = np.where(s == "+", te, ts)     # ... and to qend
flip = (dom == "-")
y_q = np.where(flip, tl - fq, fq)    # y at qstart in the oriented frame
y_qe = np.where(flip, tl - fqe, fqe)  # y at qend
df["x"] = df["qname"].map(chr_x0).values + df["qstart"].values
df["xend"] = df["qname"].map(chr_x0).values + df["qend"].values
df["y"] = df["y_off"].values + y_q
df["yend"] = df["y_off"].values + y_qe
df["syntenic"] = (df["strand"].values == dom)

if SYNTENIC_ONLY:
    kept = df["syntenic"]
    print(f"dropping {(~kept).mean():.1%} of alignments (opposite strand within a "
          f"contig) -> {int(kept.sum()):,} kept")
    df_plot = df[kept].copy()
else:
    df_plot = df.copy()

# ============================================================
# Whole-genome figure
# ============================================================
color_of = dict(zip(chr_order, plt.cm.turbo(np.linspace(0, 1, len(chr_order)))))

fig, ax = plt.subplots(figsize=(11, 9))

for c in chr_order:
    sub = df_plot[df_plot["qname"] == c]
    ax.scatter(sub["x"], sub["y"], s=DOT_SIZE, alpha=0.65, color=color_of[c],
               linewidths=0, edgecolors="none")

# chromosome boundary lines
for x0 in np.r_[0, chr_len.cumsum().values]:
    ax.axvline(x0, color="0.45", lw=0.6, ls="--", alpha=0.7)
band = {}
for t in sorted_contigs:
    c = prim_chr[t]
    if c not in band:
        band[c] = y_off[t]
for y0 in sorted(band.values()):
    ax.axhline(y0, color="0.45", lw=0.6, ls="--", alpha=0.7)

ax.set_xlim(0, GENOME_X)
ax.set_ylim(0, y_total)
xt = np.arange(0, GENOME_X, 2.5e8)
yt = np.arange(0, y_total, 2.5e8)
ax.set_xticks(xt, [f"{int(t/1e6)} Mb" for t in xt], rotation=45, fontsize=13)
ax.set_yticks(yt, [f"{int(t/1e6)} Mb" for t in yt], fontsize=13)

ax.set_xlabel(f"OctDeg2.0 assembly, concatenated ({len(chr_order)} aligned chromosomes)", fontsize=15)
ax.set_ylabel(f"OctDeg1.0 assembly, ordered & oriented ({len(sorted_contigs):,} aligned contigs)", fontsize=15)
ax.set_title(f"Whole-genome collinearity: OctDeg1.0 (female) vs OctDeg2.0\n"
             f"mapq>={MAPQ_MIN}, aln_len>={ALN_LEN_MIN:,}"
             + (", syntenic-only" if SYNTENIC_ONLY else ""),
             fontsize=17, pad=22)

handles = [mpatches.Patch(color=color_of[c], label=c) for c in chr_order]
ax.legend(handles=handles, title="OctDeg2.0\nscaffolds", ncol=1,
          fontsize=8, title_fontsize=10,
          bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0.)

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
print("saved:", OUT_PNG)
