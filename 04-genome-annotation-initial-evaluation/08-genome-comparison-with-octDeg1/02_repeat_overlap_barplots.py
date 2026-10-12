#!/usr/bin/env python3
"""
Panel B of Figure 3 — how much repetitive and non-repetitive OctDeg2.0 sequence
aligns to OctDeg1.0 (replaces the former separate panels B and C).

One bar per repeat class of OctDeg2.0, height = sequence in that class (Gb).
The dark part is the sequence aligned to OctDeg1.0, the light part is not
aligned. Two numbers are annotated, which together carry both former panels:
  * above each bar: % of that class that is aligned           (former panel C)
  * inside the dark part: % of all aligned sequence in that class (former panel B)

Computation (equivalent to figure/barplots-feature-overlaps-with-repeat/
alignment_in_repeat_region.ipynb, reimplemented in pure pandas/numpy so no
bedtools/pyranges is required):
  1. merge the aligned query-side intervals,
  2. merge the RepeatMasker intervals,
  3. exact basepair overlap via a two-pointer sweep over the merged intervals,
  4. non-repetitive genome = chromosomes with repeats minus repetitive bp.
The four basepair totals are cached in $FIG3_OUT_DIR/figure3B_repeat_alignment_bp.tsv;
delete nothing — pass --recompute to rebuild the cache (the old one is kept as *.pre-*).

Run with the python-visualizations env:
    python 02_repeat_overlap_barplots.py [--recompute] [--width-mm 58] [--height-mm 95]
"""

import argparse
import os
import subprocess
import sys
import time

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


FONT_SIZE = 8                                  # pt at final print size; Genome Research: 8-10 pt
ALIGNED = _cfg("FIG3_COLOR_OCTDEG2")           # aligned to OctDeg1.0
UNALIGNED = "#C9DDEC"                          # light tint of the same blue: not aligned

ap = argparse.ArgumentParser()
ap.add_argument("--recompute", action="store_true", help="recompute the basepair totals")
ap.add_argument("--width-mm", type=float, default=58)
ap.add_argument("--height-mm", type=float, default=95)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4

out_dir = _cfg("FIG3_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
cache = os.path.join(out_dir, "figure3B_repeat_alignment_bp.tsv")
out_stem = os.path.join(out_dir, "figure3B_repeat_alignment")


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
        as_, ae = ga["start"].values, ga["end"].values
        bs_, be = gb["start"].values, gb["end"].values
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


if opts.recompute or not os.path.exists(cache):
    aligned = pd.read_csv(_cfg("FIG3_ALIGNED_BED"), sep="\t", names=["chrom", "start", "end"])
    merged_aligned = merge_intervals(aligned)
    repeat = pd.read_csv(_cfg("FIG2B_REPEAT_GFF"), sep="\t", header=None, usecols=[0, 3, 4],
                         names=["chrom", "start", "end"])
    merged_repeat = merge_intervals(repeat)
    ov = overlap_bp(merged_aligned, merged_repeat)
    total_aligned = int((merged_aligned["end"] - merged_aligned["start"]).sum())
    total_repeat = int((merged_repeat["end"] - merged_repeat["start"]).sum())
    fai = pd.read_csv(_cfg("FIG3_CHROM_FAI"), sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"])
    total_genome = int(fai[fai["chrom"].isin(merged_repeat["chrom"].unique())]["length"].sum())
    tbl = pd.DataFrame({
        "repeat_class": ["Non-repetitive", "Repetitive"],
        "class_bp": [total_genome - total_repeat, total_repeat],
        "aligned_bp": [total_aligned - ov, ov],
    })
    if os.path.exists(cache):
        os.replace(cache, f"{cache}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    tbl.to_csv(cache, sep="\t", index=False)
tbl = pd.read_csv(cache, sep="\t")
tbl["pct_class_aligned"] = tbl["aligned_bp"] / tbl["class_bp"] * 100          # former panel C
tbl["pct_of_aligned"] = tbl["aligned_bp"] / tbl["aligned_bp"].sum() * 100     # former panel B
print(tbl.assign(class_Gb=tbl.class_bp / 1e9, aligned_Gb=tbl.aligned_bp / 1e9).round(3).to_string(index=False))

# ============================================================
# Plot
# ============================================================
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE, "ytick.labelsize": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
LEFT_IN, RIGHT_IN, TOP_IN, BOTTOM_IN = 0.40, 0.05, 0.50, 0.25
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
ax = fig.add_axes((LEFT_IN / PANEL_W, BOTTOM_IN / PANEL_H,
                   1 - (LEFT_IN + RIGHT_IN) / PANEL_W, 1 - (TOP_IN + BOTTOM_IN) / PANEL_H))
x = np.arange(len(tbl))
cls_gb, aln_gb = tbl["class_bp"] / 1e9, tbl["aligned_bp"] / 1e9
ax.bar(x, cls_gb, width=0.62, color=UNALIGNED, lw=0, zorder=2)
ax.bar(x, aln_gb, width=0.62, color=ALIGNED, lw=0, zorder=3)
for xi, c, a, pc, pa in zip(x, cls_gb, aln_gb, tbl["pct_class_aligned"], tbl["pct_of_aligned"]):
    ax.text(xi, c + 0.03, f"{pc:.1f}%\naligned", ha="center", va="bottom", linespacing=1.0)
    ax.text(xi, a / 2, f"{pa:.1f}%\nof aligned", ha="center", va="center", color="white", linespacing=1.0)
ax.set_xticks(x, tbl["repeat_class"])
ax.tick_params(axis="x", length=0, pad=3)
ax.set_xlim(-0.55, len(tbl) - 0.45)
ax.set_ylim(0, 2.0)
ax.set_yticks([0, 0.5, 1, 1.5, 2], ["0", "0.5", "1", "1.5", "2"])
ax.set_ylabel("OctDeg2.0 sequence (Gb)", labelpad=2)
ax.yaxis.grid(True, color="#E5E5E5", lw=0.5, zorder=0)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
fig.legend(handles=[Patch(color=ALIGNED, label="Aligned to OctDeg1.0"),
                    Patch(color=UNALIGNED, label="Not aligned")],
           loc="upper left", bbox_to_anchor=(LEFT_IN / PANEL_W, 1 - 0.03 / PANEL_H),
           frameon=False, ncol=1, handlelength=1.2, handleheight=0.9, borderaxespad=0,
           labelspacing=0.3)

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
