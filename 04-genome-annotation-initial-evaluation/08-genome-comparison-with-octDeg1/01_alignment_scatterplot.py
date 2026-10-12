#!/usr/bin/env python3
"""
Panel A of Figure 3 — whole-genome collinearity between the new OctDeg2.0
assembly (query; hifiasm-041425) and the OctDeg1.0 reference (target),
aligned with minimap2 (-x asm5).

Each alignment is drawn as a segment. The x-axis concatenates the OctDeg2.0
chromosomes in natural order; the y-axis orders and orients the OctDeg1.0
contigs so that collinear alignments run along the diagonal. Only
high-confidence alignments (mapq >= MAPQ_MIN, alignment length >= ALN_LEN_MIN)
are shown, and optionally only alignments on each contig's dominant strand
(syntenic).

Drawn at printed size for Genome Research (8 pt text, no title). Chromosomes
are labelled by number under the x-axis (two staggered rows) and shaded in
alternating tones of the OctDeg2.0 colour, instead of a 30-entry colour key.

Port of figure/align-octDeg1-hifiasm-minimap2/visualize_paf_clean_collinearity.ipynb.

Inputs: $FIG3_PAF (config.sh). Outputs: $FIG3_OUT_DIR/figure3A_alignment_dotplot.{pdf,svg,png}

Run with the python-visualizations env:
    python 01_alignment_scatterplot.py [--width-mm 110] [--height-mm 95]
"""

import argparse
import os
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


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


# ---------------------------------------------------------------------------
# RUN PARAMETERS
# ---------------------------------------------------------------------------
MAPQ_MIN = 60            # keep primary, high-confidence mappings only
ALN_LEN_MIN = 10_000     # bp; drop short / repetitive noise
SYNTENIC_ONLY = True     # keep only alignments matching each contig's dominant strand
FONT_SIZE = 8            # pt at final print size; Genome Research: 8-10 pt
SEG_LW = 0.8             # pt; alignment segments
COLOR_A = _cfg("FIG3_COLOR_OCTDEG2")   # OctDeg2.0 blue; alternate chromosomes use the light tint
COLOR_B = "#7FB8DE"

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=110)
ap.add_argument("--height-mm", type=float, default=95)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4

out_dir = _cfg("FIG3_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, "figure3A_alignment_dotplot")


def natural_key(s):
    """Natural sort key for chromosome names (chr1 < chr2 < ... < chr10 < chrX)."""
    return re.sub(r"(\d+)", lambda m: m.group(1).zfill(8), s)

# ============================================================
# Load + filter PAF
# ============================================================
cols = ["qname", "qlen", "qstart", "qend", "strand", "tname",
        "tlen", "tstart", "tend", "matches", "aln_len", "mapq"]
raw = pd.read_csv(_cfg("FIG3_PAF"), sep="\t", header=None, usecols=range(12), names=cols)
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
print(f"contigs drawn (y-axis): {df_plot['tname'].nunique():,}")

# ============================================================
# Figure at printed size
# ============================================================
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 0.42, 0.06, 0.05, 0.52
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
ax = fig.add_axes((LEFT_IN / PANEL_W, BOTTOM_IN / PANEL_H,
                   1 - (LEFT_IN + RIGHT_IN) / PANEL_W, 1 - (TOP_IN + BOTTOM_IN) / PANEL_H))

# chromosome boundaries (x) and OctDeg2.0-chromosome bands of OctDeg1.0 contigs (y)
for x0 in chr_len.cumsum().values[:-1]:
    ax.axvline(x0, color="#DDDDDD", lw=0.4, zorder=0)
band = {}
for t in sorted_contigs:
    band.setdefault(prim_chr[t], y_off[t])
for y0 in sorted(band.values())[1:]:
    ax.axhline(y0, color="#DDDDDD", lw=0.4, zorder=0)

for i, c in enumerate(chr_order):
    sub = df_plot[df_plot["qname"] == c]
    segs = np.stack([sub[["x", "y"]].values, sub[["xend", "yend"]].values], axis=1)
    ax.add_collection(LineCollection(segs, colors=COLOR_A if i % 2 == 0 else COLOR_B,
                                     linewidths=SEG_LW, capstyle="round", zorder=2))

ax.set_xlim(0, GENOME_X)
ax.set_ylim(0, y_total)
yt = np.arange(0, y_total, 5e8)
ax.set_yticks(yt, [f"{t/1e9:g}" for t in yt])
ax.set_ylabel("OctDeg1.0 contigs, ordered and oriented (Gb)", labelpad=2)
ax.set_xticks(chr_len.cumsum().values[:-1])
ax.tick_params(axis="x", labelbottom=False, length=2)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

# chromosome numbers under the x-axis, alternating between two rows
row_y = [-0.012, -0.012 - 0.13 / (PANEL_H - TOP_IN - BOTTOM_IN)]
for i, c in enumerate(chr_order):
    xc = (chr_x0[c] + chr_len[c] / 2) / GENOME_X
    ax.text(xc, row_y[i % 2], c.replace("chr", ""), transform=ax.transAxes,
            ha="center", va="top")
ax.text(0.5, row_y[1] - 0.11 / (PANEL_H - TOP_IN - BOTTOM_IN), "OctDeg2.0 chromosome",
        transform=ax.transAxes, ha="center", va="top")

# Within each row, push crowded labels apart (small chromosomes 24-28) to keep a ~0.8 mm gap
fig.canvas.draw()
r = fig.canvas.get_renderer()
labels = ax.texts[:len(chr_order)]
gap = 0.03 / (PANEL_W - LEFT_IN - RIGHT_IN)                      # axes fraction
for row in (labels[0::2], labels[1::2]):
    w = [t.get_window_extent(r).transformed(ax.transAxes.inverted()).width for t in row]
    xs = [t.get_position()[0] for t in row]
    for i in range(1, len(xs)):
        xs[i] = max(xs[i], xs[i - 1] + (w[i - 1] + w[i]) / 2 + gap)
    for i in range(len(xs) - 1, -1, -1):
        lim = 1.0 - w[i] / 2 if i == len(xs) - 1 else xs[i + 1] - (w[i] + w[i + 1]) / 2 - gap
        xs[i] = min(xs[i], lim)
    for t, x in zip(row, xs):
        t.set_x(x)
    assert all(b - a >= (w[i] + w[i + 1]) / 2 + gap * 0.99
               for i, (a, b) in enumerate(zip(xs, xs[1:]))), "chromosome labels overlap"

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
