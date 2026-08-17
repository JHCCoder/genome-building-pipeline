#!/usr/bin/env python
"""make_occupancy_plot_notebook.py — assemble the genome-wide occupancy plot notebook."""
import nbformat as nbf
import os

WORK = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
OUT = os.path.join(WORK, "repeat_occupancy_genomewide_plot.ipynb")

nb = nbf.v4.new_notebook()
nb["metadata"] = {
    "kernelspec": {"display_name": "Python (python-visualizations)", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.13.1"},
}

cells = []

cells.append(nbf.v4.new_markdown_cell(
"""# Genome-wide repeat-class occupancy around the CENP-A core

One connected scatterplot per repeat-period bin (small multiples, free y-scale).
For each bin the x-axis is the spatial interval ordered along the chromosome:
**short side → core → long side → background**; the y-axis is the % of that
interval's bp covered by the bin.

Because the bins span very different densities (348–349 bp reaches ~31% while
most bins are <3%), each bin gets its own y-scale — a shared axis would flatten
everything except the 349-bp satellite.

**chrY is excluded** from the genome-wide aggregate and the per-chromosome
figures: no tandem-repeat (TRF) analysis was performed on chrY, so it carries
0 repeat bp and would only dilute every percentage with a flat-zero
contribution. The genome-wide rows sum the 29 autosomes + chrX only.

**CENP-A core = MACS2 peaks.** The core is no longer a fixed 100 kb window: it is
the per-chromosome CENP-A peak called by MACS2 (CENP-A CUT&Tag reps vs H3K27ac
control, `--nomodel --shift -100 --extsize 200`, q ≤ 0.05), overlapping the
strict-signal peak anchor. Cores are typically 2–27 kb (chr4 = the 16.7 kb
degenerate junction). All other intervals are unchanged.

Input: `results/repeat_occupancy_peakbased_with_genome.csv` (genome-wide rows only).
"""))

cells.append(nbf.v4.new_code_cell(
"""import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

WORK = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
# peak-based CENP-A core (MACS2 peaks) — see 10_call_cenpa_peaks.sh / 10b
occ  = pd.read_csv(os.path.join(WORK, "results", "repeat_occupancy_peakbased_with_genome.csv"))
gw   = occ[occ["chromosome"] == "genome"].copy()

# --- ordered spatial intervals along the chromosome ---
interval_order = [
    ("short", "distal"), ("short", "intermediate"), ("short", "proximal"),
    ("core", "core"),
    ("long", "proximal"), ("long", "intermediate"), ("long", "distal"),
    ("background", "background"),
]
gw["x"] = gw.apply(lambda r: interval_order.index((r["side"], r["region"])), axis=1)
xlab = ["short\\ndistal", "short\\nintermediate", "short\\nproximal", "CENP-A\\ncore",
        "long\\nproximal", "long\\nintermediate", "long\\ndistal", "background"]

# --- the 9 TRF period bins, fixed hue order (project palette) ---
# bins 3 (51-192) and 4 (193-195) were pale cream / pale teal and nearly
# invisible on white; darkened to visible ochre and dark teal respectively.
BINS = ["1-10 bp","11-50 bp","51-192 bp","193-195 bp","196-347 bp",
        "348-349 bp","350-385 bp","386-390 bp","391+ bp"]
PAL  = {"1-10 bp":"#8c510a","11-50 bp":"#d8b365","51-192 bp":"#d9a106",
        "193-195 bp":"#35978f","196-347 bp":"#5ab4ac","348-349 bp":"#2166ac",
        "350-385 bp":"#762a83","386-390 bp":"#e7298a","391+ bp":"#1b7837"}
print("loaded", len(gw), "genome-wide rows")
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- helper: aggregate split segments per (side, region) before plotting ----
# Background (and occasionally a flank) can be split into multiple non-contiguous
# segments; sum bp per (side, region, bin) so each region is ONE point on the
# x-axis (otherwise a connected scatterplot draws a spurious vertical line and
# an extra dot, e.g. the two chr4 background segments both mapping to x=7).
def agg_region(df):
    g = (df.groupby(["side","region","region_order","repeat_bin","repeat_bin_label"], as_index=False)
           .agg(interval_bp=("interval_bp","sum"),
                repeat_overlap_bp=("repeat_overlap_bp","sum"),
                repeat_interval_count=("repeat_interval_count","sum")))
    g["repeat_percent"] = 100.0 * g["repeat_overlap_bp"] / g["interval_bp"]
    return g

# ---- small multiples: one connected scatterplot per bin (free y) ----
def draw_bin_panel(ax, d, bin_label, tick_fs=8, title_fs=10, mk=5, lw=1.6):
    d = d.sort_values("x")
    ax.plot(d["x"], d["repeat_percent"], marker="o", markersize=mk,
            linewidth=lw, color=PAL[bin_label], zorder=3)
    ax.fill_between(d["x"], d["repeat_percent"], 0, alpha=0.08, color=PAL[bin_label], zorder=1)
    for xs, xe, c in [(0,2.5,"#eaf2f8"), (3.5,6.5,"#fef6e7")]:
        ax.axvspan(xs, xe, color=c, zorder=0)
    ax.axvline(3.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)   # core | long
    ax.axvline(6.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)   # long | background
    ax.set_xticks(range(8))
    ax.set_xticklabels(xlab, fontsize=tick_fs, rotation=30, ha="right")
    ax.tick_params(axis="y", labelsize=tick_fs)
    ax.set_title(f"{bin_label}  ({d['repeat_percent'].max():.1f}% max)", fontsize=title_fs)
    ax.grid(axis="y", color="#e0e0e0", lw=0.5, zorder=0)
    ax.spines[["top","right"]].set_visible(False)

gw_a = agg_region(gw)
gw_a["x"] = gw_a.apply(lambda r: interval_order.index((r["side"], r["region"])), axis=1)

# ---- general small-multiples figure from any subset of bins ----
def draw_occupancy_figure(bin_list, title, fname, ncols=3, annotate_sides=True,
                          tick_fs=8, title_fs=10, side_fs=7.5, suptitle_fs=13,
                          panel_w=6.5, panel_h=4.6, mk=5, lw=1.6):
    n = len(bin_list)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_w*ncols, panel_h*nrows))
    axes = axes.ravel()
    for ax, bin_label in zip(axes, bin_list):
        draw_bin_panel(ax, gw_a[gw_a["repeat_bin_label"] == bin_label], bin_label,
                       tick_fs=tick_fs, title_fs=title_fs, mk=mk, lw=lw)
    for ax in axes[n:]:
        ax.set_visible(False)
    if annotate_sides:
        for ax in axes[:n]:
            ytop = ax.get_ylim()[1]
            ax.text(1.25, ytop*0.97, "SHORT SIDE", ha="center", va="top",
                    fontsize=side_fs, color="#4a7ba6", fontweight="bold")
            ax.text(5, ytop*0.97, "LONG SIDE", ha="center", va="top",
                    fontsize=side_fs, color="#a67c2e", fontweight="bold")
    fig.suptitle(title, fontsize=suptitle_fs)
    fig.tight_layout(rect=[0,0,1,0.96])
    outdir = os.path.join(WORK, "plots")
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, fname + ".png")
    fig.savefig(png, dpi=150)
    fig.savefig(os.path.join(outdir, fname + ".pdf"))
    plt.close(fig)
    from IPython.display import Image, display
    display(Image(filename=png))

# 1) all nine bins
draw_occupancy_figure(BINS,
    "Genome-wide repeat-class occupancy by distance from the CENP-A core\\n"
    "(% of interval bp; one panel per TRF period bin, free y-scale)",
    "genomewide_occupancy_by_bin", ncols=3)
print("saved plots/genomewide_occupancy_by_bin.{png,pdf}")

# 2) key bins: 1 (microsat), 4 (195 = L1 5' tandem), 6 (349), 8 (389) — LARGE FONT
SUBSET_A = ["1-10 bp", "193-195 bp", "348-349 bp", "386-390 bp"]
draw_occupancy_figure(SUBSET_A,
    "Genome-wide occupancy, key bins: 1-10 bp, 193-195 bp, 348-349 bp, 386-390 bp\\n"
    "(% of interval bp; free y-scale)",
    "genomewide_occupancy_key_bins", ncols=2,
    tick_fs=15, title_fs=19, side_fs=14, suptitle_fs=20,
    panel_w=9.5, panel_h=6.5, mk=8, lw=2.4)

# 3) the rest — LARGE FONT
SUBSET_B = [b for b in BINS if b not in SUBSET_A]
draw_occupancy_figure(SUBSET_B,
    "Genome-wide occupancy, remaining bins: 11-50, 51-192, 196-347, 350-385, 391+ bp\\n"
    "(% of interval bp; free y-scale)",
    "genomewide_occupancy_rest_bins", ncols=3,
    tick_fs=14, title_fs=18, side_fs=13, suptitle_fs=19,
    panel_w=8.5, panel_h=5.8, mk=7, lw=2.2)
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- legend: the 9 bins in fixed hue order ----
import io as _io
legend_handles = [Line2D([0],[0], marker="o", ls="-", color=PAL[b], label=b) for b in BINS]
fig, ax = plt.subplots(figsize=(9, 0.8))
ax.legend(handles=legend_handles, loc="center", ncol=9, frameon=False, fontsize=10)
ax.axis("off")
_buf = _io.BytesIO()
fig.savefig(_buf, format="png", dpi=150)
_buf.seek(0)
from IPython.display import Image as _LImg, display as _Ldisp
_Ldisp(_LImg(data=_buf.getvalue()))
plt.close(fig)
"""))

cells.append(nbf.v4.new_markdown_cell(
"""### Per-chromosome figures

The same 9-bin small-multiples layout, produced per chromosome with a plotting
loop. All 30 chromosomes are written to `plots/per_chromosome/occupancy_<chr>.pdf`
and `.png`; a few representative ones are shown inline below (chr4 = clean 349-bp
array spanning the short side + core; chr1 = 349 in short-side flanks but not core;
chr9 = 349 through core + flanks; chrX = no 349 near the peak).
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- per-chromosome plotting loop (displays inline in the notebook) ----
# One small-multiples figure per chromosome (same 9-bin layout as genome-wide,
# free y per panel). Figures are shown inline in the notebook; PDF/PNG are also
# written to plots/per_chromosome/ for the paper.
import os, io
from IPython.display import Image as _Image, display as _display

def plot_one_chromosome(chrom):
    d = occ[occ["chromosome"] == chrom].copy()
    if d.empty:
        return None
    d = agg_region(d)
    d["x"] = d.apply(lambda r: interval_order.index((r["side"], r["region"])), axis=1)
    fig, axes = plt.subplots(3, 3, figsize=(15, 11))
    axes = axes.ravel()
    for ax, bin_label in zip(axes, BINS):
        dd = d[d["repeat_bin_label"] == bin_label]
        if dd.empty:
            ax.set_visible(False); continue
        draw_bin_panel(ax, dd, bin_label)
    for ax in axes[3:6]:
        if not ax.get_visible(): continue
        ax.text(1.25, ax.get_ylim()[1]*0.97, "SHORT SIDE", ha="center", va="top",
                fontsize=8, color="#4a7ba6", fontweight="bold")
        ax.text(5, ax.get_ylim()[1]*0.97, "LONG SIDE", ha="center", va="top",
                fontsize=8, color="#a67c2e", fontweight="bold")
    fig.suptitle(f"{chrom}  repeat-class occupancy by distance from CENP-A core\\n"
                 f"(% of interval bp; one panel per TRF period bin, free y-scale)", fontsize=13)
    fig.tight_layout(rect=[0,0,1,0.96])
    outdir = os.path.join(WORK, "plots", "per_chromosome")
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"occupancy_{chrom}.png")
    fig.savefig(png, dpi=110)
    fig.savefig(os.path.join(outdir, f"occupancy_{chrom}.pdf"))
    # render to buffer and display inline (works with Agg backend)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    buf.seek(0)
    _display(_Image(data=buf.getvalue()))
    plt.close(fig)
    return png

# chrY is excluded: no tandem-repeat (TRF) analysis was performed on chrY, so
# every repeat_overlap_bp would be 0 and its figures would be flat-zero noise.
per_chr_chroms = [c for c in occ["chromosome"].unique() if c not in ("genome", "chrY")]
print(f"plotting {len(per_chr_chroms)} chromosomes (chrY excluded; inline + plots/per_chromosome/)")
for chrom in per_chr_chroms:
    plot_one_chromosome(chrom)
print("done")
"""))

cells.append(nbf.v4.new_markdown_cell(
"""### HiFi read coverage around the CENP-A core

The same spatial-interval philosophy applied to **HiFi (long-read) coverage**.
Source: `contig-coverage/2-primaryreads-coverage/..._chrAssigned_primaryReads_15kb_windows.tsv`
— mean HiFi coverage per 15 kb window, **primary alignments only** (i.e. uniquely
mapped), on the chr-assigned assembly whose chromosome lengths match this
analysis (chr1 = 185,300,982). chrY is excluded to match the repeat analysis.
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- shared setup for the HiFi section: chromosome sizes, domains, intervals ----
import io as _io
from IPython.display import Image as _Image, display as _display
SRC_C = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
CORE_WIN = 100_000

sizes = {}
with open(os.path.join(SRC_C, "data", "chrom_sizes.txt")) as _f:
    for _line in _f:
        _p = _line.split()
        if _p[0].startswith("chr"):
            sizes[_p[0]] = int(_p[1])
CHROMS = sorted(sizes, key=lambda c: (c[3:].isdigit(), int(c[3:]) if c[3:].isdigit() else (1 if c=="X" else 2)))

# PEAK-BASED CENP-A cores (MACS2) — see 10_call_cenpa_peaks.sh / 10b_peaks_to_cores.py
dom = pd.read_csv(os.path.join(WORK, "data", "domains", "cenpa_domains_peakbased.csv"))
dom = dom[dom["chrom"].isin(CHROMS)].copy()
dom["chrom_len"] = dom["chrom"].map(sizes)
dom["core_start"] = dom["core_start"].astype(int)
dom["core_end"] = np.minimum(dom["core_end"].astype(int), dom["chrom_len"])

# re-create interval geometry helpers used by the occupancy frame
EDGES = [("proximal", 0, 500_000),
         ("intermediate", 500_000, 2_000_000),
         ("distal", 2_000_000, 5_000_000)]
REGION_ORDER = {"short_distal":1, "short_intermediate":2, "short_proximal":3,
                "core":4, "long_proximal":5, "long_intermediate":6, "long_distal":7,
                "background":8}

def _fl_left(cs, lo, hi):
    s = max(0, cs - hi); e = max(0, cs - lo)
    return (s, e) if e > s else None

def _fl_right(ce, clen, lo, hi):
    s = min(clen, ce + lo); e = min(clen, ce + hi)
    return (s, e) if e > s else None

def build_intervals(cs, ce, clen):
    left_len, right_len = cs, clen - ce
    short_side = "left" if left_len <= right_len else "right"
    long_side = "right" if short_side == "left" else "left"
    ivs = []
    for (region, lo, hi) in EDGES:
        for side in ("short", "long"):
            arm = short_side if side == "short" else long_side
            iv = _fl_left(cs, lo, hi) if arm == "left" else _fl_right(ce, clen, lo, hi)
            if iv: ivs.append((side, region, *iv))
    ivs.append(("core", "core", cs, ce))
    excl = sorted((s, e) for (_, _, s, e) in ivs)
    bg, cur = [], 0
    for (s, e) in excl:
        if s > cur: bg.append((cur, s))
        cur = max(cur, e)
    if cur < clen: bg.append((cur, clen))
    for (s, e) in bg:
        if e > s: ivs.append(("background", "background", s, e))
    return ivs, short_side, long_side

interval_order = [
    ("short","distal"), ("short","intermediate"), ("short","proximal"),
    ("core","core"), ("long","proximal"), ("long","intermediate"),
    ("long","distal"), ("background","background")]
print("setup ready:", len(CHROMS), "chromosomes,", len(dom), "domains")
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- load HiFi primary-read coverage (15 kb windows, chr1-28,X,Y) ----
HIFI_TSV = ("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/"
            "code/command-line-script/contig-coverage/2-primaryreads-coverage/"
            "coverage_long_read_hifiasm_041425_scaffolded_juiceBox_sorted_"
            "chrAssigned_primaryReads_15kb_windows.tsv")
hifi = pd.read_csv(HIFI_TSV, sep="\\t", header=None,
                   names=["chrom", "start", "end", "mean_cov"])
hifi = hifi[hifi["chrom"].isin(CHROMS)].copy()
hifi["mid"] = (hifi["start"] + hifi["end"]) / 2
print(f"HiFi windows (chr1-28,X,Y): {len(hifi):,}; genome mean coverage: {hifi['mean_cov'].mean():.1f}x")
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- HiFi coverage per (chromosome, spatial interval): MEAN and MAX ----
# Mean captures the bulk (array-average) coverage; MAX captures the collapse
# hotspot (reads piling at the one uniquely-mappable locus in a collapsed array).
dom_c = dom[["chrom", "peak_window", "core_start", "core_end"]].copy()
hifi_rows = []
for _, d in dom_c.iterrows():
    chrom = d["chrom"]; cs = int(d["core_start"]); ce = int(d["core_end"])
    clen = int(sizes[chrom])
    ivs, _, _ = build_intervals(cs, ce, clen)
    sub = hifi[hifi["chrom"] == chrom]
    for (side, region, s, e) in ivs:
        ibp = e - s
        if ibp <= 0:
            continue
        win = sub[(sub["start"] < e) & (sub["end"] > s)].copy()
        if len(win) == 0:
            continue
        win["ov"] = np.maximum(0, np.minimum(win["end"], e) - np.maximum(win["start"], s))
        mean_cov = (win["mean_cov"] * win["ov"]).sum() / win["ov"].sum()
        max_cov  = win["mean_cov"].max()
        region_key = f"{side}_{region}" if side not in ("core","background") else side
        hifi_rows.append({"chromosome": chrom, "side": side, "region": region,
                          "region_order": REGION_ORDER[region_key],
                          "interval_start": s, "interval_end": e, "interval_bp": ibp,
                          "mean_hifi_cov": mean_cov, "max_hifi_cov": max_cov})
hifi_occ = pd.DataFrame(hifi_rows)
print(f"per-chromosome HiFi interval rows: {len(hifi_occ)}")
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- genome-wide HiFi coverage per interval (bp-summed across chr, chrY excl) ----
# MAX coverage per interval is the collapse signature: within a collapsed array,
# reads pile at the single uniquely-mappable locus, giving a coverage hotspot
# (e.g. chr4 349-array junction 154.7x vs chr mean 35x). MEAN is shown as a
# reference (bulk array coverage, often below background due to ambiguous mapping).
agg_f = lambda g: pd.Series({
    "mean_hifi_cov": (g["mean_hifi_cov"] * g["interval_bp"]).sum() / g["interval_bp"].sum(),
    "max_hifi_cov": g["max_hifi_cov"].max(),
    "interval_bp": g["interval_bp"].sum()})
hifi_gw = (hifi_occ[hifi_occ["chromosome"] != "chrY"]
    .groupby(["side", "region", "region_order"], as_index=False).apply(agg_f, include_groups=False))
hifi_gw["x"] = hifi_gw.apply(lambda r: interval_order.index((r["side"], r["region"])), axis=1)
hifi_gw = hifi_gw.sort_values("x")

# ---- plot: genome-wide HiFi coverage by interval (connected scatterplot) ----
fig, ax = plt.subplots(figsize=(10, 7.0))
ax.plot(hifi_gw["x"], hifi_gw["max_hifi_cov"], marker="o", markersize=8,
        linewidth=2.2, color="#8c3a1f", zorder=3, label="max coverage (collapse hotspot)")
ax.plot(hifi_gw["x"], hifi_gw["mean_hifi_cov"], marker="s", markersize=6,
        linewidth=1.6, color="#2166ac", zorder=3, label="mean coverage (bulk)")
ax.fill_between(hifi_gw["x"], hifi_gw["mean_hifi_cov"], 0, alpha=0.08, color="#2166ac", zorder=1)
for xs, xe, c in [(0,2.5,"#eaf2f8"), (3.5,6.5,"#fef6e7")]:
    ax.axvspan(xs, xe, color=c, zorder=0)
ax.axvline(3.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)
ax.axvline(6.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)
ax.set_xticks(range(8))
ax.set_xticklabels(xlab, fontsize=12, rotation=30, ha="right")
ax.tick_params(axis="y", labelsize=12)
ax.set_title("Genome-wide HiFi (primary-read) coverage by distance from CENP-A core\\n"
             "(max & mean coverage per interval; chrY excluded)", fontsize=14)
ax.set_ylabel("HiFi coverage (x)", fontsize=13)
ax.grid(axis="y", color="#e0e0e0", lw=0.5, zorder=0)
ax.spines[["top","right"]].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 0.945),
          bbox_transform=fig.transFigure, ncol=2, fontsize=11, frameon=False)
yt = ax.get_ylim()[1]
ax.text(1.25, yt*0.97, "SHORT SIDE", ha="center", va="top", fontsize=11,
        color="#4a7ba6", fontweight="bold")
ax.text(5, yt*0.97, "LONG SIDE", ha="center", va="top", fontsize=11,
        color="#a67c2e", fontweight="bold")
fig.subplots_adjust(left=0.08, right=0.96, top=0.87, bottom=0.11)
_buf = io.BytesIO(); fig.savefig(_buf, format="png", dpi=150); _buf.seek(0)
_display(_Image(data=_buf.getvalue()))
fig.savefig(os.path.join(WORK, "plots", "genomewide_hifi_coverage_by_interval.png"), dpi=150)
fig.savefig(os.path.join(WORK, "plots", "genomewide_hifi_coverage_by_interval.pdf"))
plt.close(fig)
print(hifi_gw[["side","region","mean_hifi_cov","max_hifi_cov"]].round(2).to_string(index=False))
"""))

cells.append(nbf.v4.new_markdown_cell(
"""### Secondary QC: all-reads HiFi coverage (unfiltered)

The primary-reads coverage above keeps only each read's primary alignment
(`-F 2304` → drop secondary + supplementary). As a QC check, the same intervals
are also evaluated with the **all-reads** alignment (no such filter): every
placement is counted, so in repetitive/low-mappability sequence the coverage is
inflated by multimapping reads (e.g. the chr4 349-array junction reads 5,504×
all-reads vs 154.7× primary). Comparing the two shows where the primary-only
filter changes the picture.

Source: `contig-coverage/1-allreads-coverage/..._chrAssigned_15kb_windows.tsv`.
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- load ALL-READS HiFi coverage (unfiltered, 15 kb windows) ----
ALL_TSV = ("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/"
           "code/command-line-script/contig-coverage/1-allreads-coverage/"
           "coverage_long_read_hifiasm_041425_scaffolded_juiceBox_sorted_"
           "chrAssigned_15kb_windows.tsv")
hifi_all = pd.read_csv(ALL_TSV, sep="\\t", header=None,
                       names=["chrom", "start", "end", "mean_cov"])
hifi_all = hifi_all[hifi_all["chrom"].isin(CHROMS)].copy()
print(f"all-reads windows (chr1-28,X,Y): {len(hifi_all):,}; "
      f"genome mean: {hifi_all['mean_cov'].mean():.1f}x")
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- all-reads coverage per (chromosome, spatial interval): mean & max ----
all_rows = []
for _, d in dom.iterrows():
    chrom = d["chrom"]; cs = int(d["core_start"]); ce = int(d["core_end"])
    clen = int(d["chrom_len"])
    ivs, _, _ = build_intervals(cs, ce, clen)
    sub = hifi_all[hifi_all["chrom"] == chrom]
    for (side, region, s, e) in ivs:
        ibp = e - s
        if ibp <= 0:
            continue
        win = sub[(sub["start"] < e) & (sub["end"] > s)].copy()
        if len(win) == 0:
            continue
        win["ov"] = np.maximum(0, np.minimum(win["end"], e) - np.maximum(win["start"], s))
        region_key = f"{side}_{region}" if side not in ("core","background") else side
        all_rows.append({"chromosome": chrom, "side": side, "region": region,
                         "region_order": REGION_ORDER[region_key],
                         "interval_start": s, "interval_end": e, "interval_bp": ibp,
                         "mean_cov": (win["mean_cov"] * win["ov"]).sum() / win["ov"].sum(),
                         "max_cov": win["mean_cov"].max()})
all_occ = pd.DataFrame(all_rows)

# ---- genome-wide all-reads coverage per interval (bp-summed, chrY excl) ----
all_gw = (all_occ[all_occ["chromosome"] != "chrY"]
    .groupby(["side", "region", "region_order"], as_index=False)
    .apply(lambda g: pd.Series({
        "mean_cov": (g["mean_cov"] * g["interval_bp"]).sum() / g["interval_bp"].sum(),
        "max_cov": g["max_cov"].max(), "interval_bp": g["interval_bp"].sum()}), include_groups=False))
all_gw["x"] = all_gw.apply(lambda r: interval_order.index((r["side"], r["region"])), axis=1)
all_gw = all_gw.sort_values("x")

# ---- plot: all-reads vs primary-reads mean coverage per interval ----
fig, ax = plt.subplots(figsize=(10, 7.0))
ax.plot(all_gw["x"], all_gw["mean_cov"], marker="o", markersize=8, linewidth=2.2,
        color="#8a5a00", zorder=3, label="all-reads (unfiltered) mean")
ax.plot(hifi_gw["x"], hifi_gw["mean_hifi_cov"], marker="s", markersize=7, linewidth=1.8,
        color="#2166ac", zorder=3, label="primary-reads mean")
for xs, xe, c in [(0,2.5,"#eaf2f8"), (3.5,6.5,"#fef6e7")]:
    ax.axvspan(xs, xe, color=c, zorder=0)
ax.axvline(3.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)
ax.axvline(6.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)
ax.set_xticks(range(8))
ax.set_xticklabels(xlab, fontsize=12, rotation=30, ha="right")
ax.tick_params(axis="y", labelsize=12)
ax.set_title("Secondary QC: all-reads vs primary-reads HiFi coverage by interval\\n"
             "(mean coverage per interval; chrY excluded)", fontsize=14)
ax.set_ylabel("mean HiFi coverage (x)", fontsize=13)
ax.grid(axis="y", color="#e0e0e0", lw=0.5, zorder=0)
ax.spines[["top","right"]].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 0.945),
          bbox_transform=fig.transFigure, ncol=2, fontsize=11, frameon=False)
yt = ax.get_ylim()[1]
ax.text(1.25, yt*0.97, "SHORT SIDE", ha="center", va="top", fontsize=11,
        color="#4a7ba6", fontweight="bold")
ax.text(5, yt*0.97, "LONG SIDE", ha="center", va="top", fontsize=11,
        color="#a67c2e", fontweight="bold")
fig.subplots_adjust(left=0.08, right=0.96, top=0.87, bottom=0.11)
_buf = io.BytesIO(); fig.savefig(_buf, format="png", dpi=150); _buf.seek(0)
_display(_Image(data=_buf.getvalue()))
fig.savefig(os.path.join(WORK, "plots", "genomewide_hifi_allreads_vs_primary.png"), dpi=150)
fig.savefig(os.path.join(WORK, "plots", "genomewide_hifi_allreads_vs_primary.pdf"))
plt.close(fig)
print("all-reads vs primary (mean coverage per interval):")
print(pd.DataFrame({"side": all_gw["side"], "region": all_gw["region"],
                    "all_reads_mean": all_gw["mean_cov"].round(2),
                    "primary_mean": hifi_gw.set_index(["side","region"])["mean_hifi_cov"]
                                    .reindex(list(zip(all_gw["side"], all_gw["region"]))).round(2).values}).to_string(index=False))
"""))

cells.append(nbf.v4.new_markdown_cell(
"""### CENP-A coverage around the CENP-A core

The same peak-based intervals, evaluated for **CENP-A signal** (mean of the two
CENP-A CUT&Tag replicates, 1/NH-weighted k=100 coverage averaged per interval).
This is the direct functional readout of CENP-A binding at and around the core.
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- CENP-A coverage per (chromosome, spatial interval) ----
# Use the 1/NH-weighted k=100 window coverage (mean of XG_150/151), same
# interval geometry as the HiFi panel (peak-based cores).
cen_w150 = pd.read_csv(os.path.join(WORK, "data", "coverage_weighted", "XG_150_win100kb_weighted.tsv"),
                       sep="\\t", header=0, names=["chrom","start","end","s150"])
cen_w151 = pd.read_csv(os.path.join(WORK, "data", "coverage_weighted", "XG_151_win100kb_weighted.tsv"),
                       sep="\\t", header=0, names=["chrom","start","end","s151"])
cen = cen_w150.merge(cen_w151, on=["chrom","start","end"])
cen["cenpa"] = (cen["s150"] + cen["s151"]) / 2
cen = cen[cen["chrom"].isin(CHROMS)].copy()

cen_rows = []
for _, d in dom.iterrows():
    chrom = d["chrom"]; cs = int(d["core_start"]); ce = int(d["core_end"])
    clen = int(d["chrom_len"])
    ivs, _, _ = build_intervals(cs, ce, clen)
    sub = cen[cen["chrom"] == chrom]
    for (side, region, s, e) in ivs:
        ibp = e - s
        if ibp <= 0:
            continue
        win = sub[(sub["start"] < e) & (sub["end"] > s)].copy()
        if len(win) == 0:
            continue
        win["ov"] = np.maximum(0, np.minimum(win["end"], e) - np.maximum(win["start"], s))
        mean_cenpa = (win["cenpa"] * win["ov"]).sum() / win["ov"].sum()
        region_key = f"{side}_{region}" if side not in ("core","background") else side
        cen_rows.append({"chromosome": chrom, "side": side, "region": region,
                         "region_order": REGION_ORDER[region_key],
                         "interval_start": s, "interval_end": e, "interval_bp": ibp,
                         "mean_cenpa_cov": mean_cenpa})
cen_occ = pd.DataFrame(cen_rows)

# ---- genome-wide mean CENP-A coverage per interval (bp-summed across chr, chrY excl) ----
cen_gw = (cen_occ[cen_occ["chromosome"] != "chrY"]
    .groupby(["side", "region", "region_order"], as_index=False)
    .apply(lambda g: pd.Series({
        "mean_cenpa_cov": (g["mean_cenpa_cov"] * g["interval_bp"]).sum() / g["interval_bp"].sum(),
        "interval_bp": g["interval_bp"].sum()}), include_groups=False))
cen_gw["x"] = cen_gw.apply(lambda r: interval_order.index((r["side"], r["region"])), axis=1)
cen_gw = cen_gw.sort_values("x")

# ---- plot: genome-wide CENP-A coverage by interval (connected scatterplot) ----
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(cen_gw["x"], cen_gw["mean_cenpa_cov"], marker="o", markersize=8,
        linewidth=2.2, color="#2166ac", zorder=3)
ax.fill_between(cen_gw["x"], cen_gw["mean_cenpa_cov"], 0, alpha=0.10, color="#2166ac", zorder=1)
for xs, xe, c in [(0,2.5,"#eaf2f8"), (3.5,6.5,"#fef6e7")]:
    ax.axvspan(xs, xe, color=c, zorder=0)
ax.axvline(3.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)
ax.axvline(6.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)
ax.set_xticks(range(8))
ax.set_xticklabels(xlab, fontsize=12, rotation=30, ha="right")
ax.tick_params(axis="y", labelsize=12)
ax.set_title("Genome-wide CENP-A coverage by distance from CENP-A core\\n"
             "(1/NH-weighted k=100, mean of reps, per interval; chrY excluded)", fontsize=14)
ax.set_ylabel("mean CENP-A coverage (weighted)", fontsize=13)
ax.grid(axis="y", color="#e0e0e0", lw=0.5, zorder=0)
ax.spines[["top","right"]].set_visible(False)
yt = ax.get_ylim()[1]
ax.text(1.25, yt*0.97, "SHORT SIDE", ha="center", va="top", fontsize=11,
        color="#4a7ba6", fontweight="bold")
ax.text(5, yt*0.97, "LONG SIDE", ha="center", va="top", fontsize=11,
        color="#a67c2e", fontweight="bold")
fig.tight_layout()
_buf = io.BytesIO(); fig.savefig(_buf, format="png", dpi=150); _buf.seek(0)
_display(_Image(data=_buf.getvalue()))
fig.savefig(os.path.join(WORK, "plots", "genomewide_cenpa_coverage_by_interval.png"), dpi=150)
fig.savefig(os.path.join(WORK, "plots", "genomewide_cenpa_coverage_by_interval.pdf"))
plt.close(fig)
print(cen_gw[["side","region","mean_cenpa_cov"]].round(2).to_string(index=False))
"""))

cells.append(nbf.v4.new_markdown_cell(
"""### Per-chromosome coverage plots

For each chromosome, three coverage tracks plotted over the same peak-based
spatial intervals (short distal → core → long distal → background): **HiFi
primary-read coverage** (mean, blue), **HiFi raw/all-reads coverage** (mean,
green; unfiltered, so multimapping inflation at collapsed arrays is visible),
and **CENP-A coverage** (1/NH-weighted k=100 mean of reps, orange). Each panel
is genome-relative but per-chromosome, so the core gradient is visible
chromosome by chromosome.
"""))

cells.append(nbf.v4.new_code_cell(
"""# ---- per-chromosome coverage plot: CENP-A + HiFi primary + HiFi raw ----
def plot_chrom_coverage(chrom):
    cs = int(dom[dom["chrom"] == chrom]["core_start"].iloc[0])
    ce = int(dom[dom["chrom"] == chrom]["core_end"].iloc[0])
    clen = int(sizes[chrom])
    ivs, short_side, long_side = build_intervals(cs, ce, clen)
    hsub = hifi[hifi["chrom"] == chrom]
    asub = hifi_all[hifi_all["chrom"] == chrom]
    csub = cen[cen["chrom"] == chrom]
    # Aggregate per (side, region): a region (notably background, sometimes a
    # flank) can split into multiple non-contiguous segments; combine them into
    # ONE point with bp-weighted means so there are no duplicate x positions.
    agg = {}
    for (side, region, s, e) in ivs:
        ibp = e - s
        if ibp <= 0:
            continue
        region_key = f"{side}_{region}" if side not in ("core","background") else side
        hv = hsub[(hsub["start"] < e) & (hsub["end"] > s)]
        av = asub[(asub["start"] < e) & (asub["end"] > s)]
        cv = csub[(csub["start"] < e) & (csub["end"] > s)]
        rec = agg.setdefault(region_key, {"x": interval_order.index((side, region)),
                                          "h_bp": 0.0, "h_sx": 0.0,
                                          "a_bp": 0.0, "a_sx": 0.0,
                                          "c_bp": 0.0, "c_sx": 0.0})
        if len(hv) > 0:
            hv = hv.copy(); hv["ov"] = np.minimum(hv["end"], e) - np.maximum(hv["start"], s)
            rec["h_bp"] += hv["ov"].sum(); rec["h_sx"] += (hv["mean_cov"] * hv["ov"]).sum()
        if len(av) > 0:
            av = av.copy(); av["ov"] = np.minimum(av["end"], e) - np.maximum(av["start"], s)
            rec["a_bp"] += av["ov"].sum(); rec["a_sx"] += (av["mean_cov"] * av["ov"]).sum()
        if len(cv) > 0:
            cv = cv.copy(); cv["ov"] = np.minimum(cv["end"], e) - np.maximum(cv["start"], s)
            rec["c_bp"] += cv["ov"].sum(); rec["c_sx"] += (cv["cenpa"] * cv["ov"]).sum()
    pts = pd.DataFrame([{
        "x": rec["x"], "side": k.split("_")[0] if "_" in k else k,
        "region": k, "interval_bp": max(rec["h_bp"], rec["c_bp"], rec["a_bp"]),
        "hifi": (rec["h_sx"] / rec["h_bp"]) if rec["h_bp"] > 0 else np.nan,
        "raw":  (rec["a_sx"] / rec["a_bp"]) if rec["a_bp"] > 0 else np.nan,
        "cenpa": (rec["c_sx"] / rec["c_bp"]) if rec["c_bp"] > 0 else np.nan,
    } for k, rec in agg.items()]).sort_values("x")

    fig, ax1 = plt.subplots(figsize=(10, 4.6))
    ax1.plot(pts["x"], pts["cenpa"], marker="o", markersize=6, linewidth=1.8,
             color="#b4542a", zorder=3, label="CENP-A coverage")
    ax1.set_ylabel("CENP-A coverage (weighted)", fontsize=12, color="#b4542a")
    ax1.tick_params(axis="y", labelcolor="#b4542a")
    ax1.set_ylim(bottom=0)   # start at zero
    # Both right axes are created as siblings from ax1 (NOT ax2 = ax1.twinx()
    # followed by ax3 = ax2.twinx()): a chained twinx misplaces the first
    # axis's tick labels to the LEFT of the plot, away from its spine.
    ax2 = ax1.twinx()
    ax3 = ax1.twinx()
    ax2.plot(pts["x"], pts["hifi"], marker="s", markersize=6, linewidth=1.8,
             color="#2166ac", zorder=3, label="HiFi primary coverage")
    ax2.set_ylabel("HiFi primary coverage (x)", fontsize=12, color="#2166ac")
    ax2.tick_params(axis="y", direction="out", labelcolor="#2166ac")
    ax2.set_ylim(bottom=0)   # start at zero
    # raw/all-reads HiFi coverage goes on the INNER right axis; at collapsed
    # arrays (e.g. chr4) it reaches ~3,500x vs ~85x primary, so it cannot share
    # ax2's scale. Primary is kept on the OUTERMOST right axis.
    ax2.spines["right"].set_position(("outward", 56))
    ax3.plot(pts["x"], pts["raw"], marker="^", markersize=6, linewidth=1.8,
             color="#1b7837", zorder=3, label="HiFi raw coverage")
    ax3.set_ylabel("HiFi raw coverage (x)", fontsize=12, color="#1b7837")
    ax3.tick_params(axis="y", direction="out", labelcolor="#1b7837")
    ax3.set_ylim(bottom=0)   # start at zero
    for xs, xe, c in [(0,2.5,"#eaf2f8"), (3.5,6.5,"#fef6e7")]:
        ax1.axvspan(xs, xe, color=c, zorder=0)
    ax1.axvline(3.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)
    ax1.axvline(6.5, color="#8a8a8a", lw=0.8, ls=":", zorder=2)
    ax1.set_xticks(range(8))
    ax1.set_xticklabels(xlab, fontsize=9, rotation=30, ha="right")
    ax1.set_title(f"{chrom}:{cs:,}-{ce:,}", fontsize=12)
    ax1.grid(axis="y", color="#e0e0e0", lw=0.5, zorder=0)
    ax1.spines[["top"]].set_visible(False)
    ax2.spines[["top"]].set_visible(False)
    ax3.spines[["top"]].set_visible(False)
    l1, lb1 = ax1.get_legend_handles_labels()
    l2, lb2 = ax2.get_legend_handles_labels()
    l3, lb3 = ax3.get_legend_handles_labels()
    ax1.legend(l1 + l2 + l3, lb1 + lb2 + lb3, loc="lower center",
               bbox_to_anchor=(0.5, 1.12), ncol=3, fontsize=10, frameon=False)
    fig.subplots_adjust(left=0.09, right=0.78, top=0.84, bottom=0.17)
    _buf = io.BytesIO(); fig.savefig(_buf, format="png", dpi=130); _buf.seek(0)
    _display(_Image(data=_buf.getvalue()))
    outdir = os.path.join(WORK, "plots", "per_chromosome")
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, f"coverage_{chrom}.png"), dpi=130)
    fig.savefig(os.path.join(outdir, f"coverage_{chrom}.pdf"))
    plt.close(fig)

per_chr_cov = [c for c in CHROMS if c != "chrY"]
print(f"plotting per-chromosome coverage for {len(per_chr_cov)} chromosomes")
for chrom in per_chr_cov:
    plot_chrom_coverage(chrom)
print("done")
"""))

cells.append(nbf.v4.new_markdown_cell(
"""### Reading the figures

* The **x-axis** walks along the chromosome: the three short-side intervals
  (distal → proximal, i.e. toward the core), the **CENP-A core**, the three
  long-side intervals (proximal → distal), then **background**.
* The **348–349 bp** panel shows the headline: ~29–33% of the short-side
  flanking bp is 349-bp satellite vs ~2.5% of background — the satellite is
  strongly oriented toward the near-end side of the CENP-A core.
* The **1–10 bp** panel shows microsatellites peaking in the core (~14%).
* The second and third figures are the same layout split into a "key bins"
  panel (1–10, 193–195, 348–349, 386–390 bp) and a "remaining bins" panel
  (11–50, 51–192, 196–347, 350–385, 391+ bp), so the low-density bins can be
  read without sharing a row with the 349-bp satellite.
* Every other bin is essentially flat (<3%), consistent with the 195-bp /
  389-bp being L1-associated and dispersed.

### Reading the HiFi coverage figure

Two lines per interval, both genome-wide (bp-summed across chr1-28,X; chrY excluded):

* **Max coverage** (orange, squares) = the collapse hotspot. Where a centromeric
  array collapsed during assembly, reads pile at the single uniquely-mappable
  locus, so the *maximum* coverage within the interval spikes. The elevated max
  in the CENP-A core (e.g. chr4 349-array junction 154.7× vs chr mean 35×) is
  the assembly-collapse signature: the satellite is under-assembled, and the
  core peaks on whatever uniquely-annealable sequence assembled.
* **Mean coverage** (blue, circles) = the bulk array coverage, often *below*
  background because reads within a homogeneous satellite map ambiguously and
  spread across copies. So mean ↓ at the core is not a contradiction of collapse
  — it is the expected dual signature (bulk low, hotspot high).

This mirrors the prior conclusion: chr4's 349 array "assembled at ~35× HiFi"
(bulk), yet its CENP-A peak sits at a ~155× collapse hotspot; on most other
chromosomes the array collapsed outright, and CENP-A piles on flanking sequence.
"""))

nb["cells"] = cells
with open(OUT, "w") as f:
    nbf.write(nb, f)
print("wrote", OUT)
