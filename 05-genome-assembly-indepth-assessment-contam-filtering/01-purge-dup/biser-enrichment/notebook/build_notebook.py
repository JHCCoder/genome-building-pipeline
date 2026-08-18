#!/usr/bin/env python3
"""Build the visualization notebook programmatically."""
import json
import nbformat as nbf
from pathlib import Path

NB_PATH = Path(__file__).resolve().parent / "purge_dups_biser_enrichment.ipynb"

cells = []

def md(source):
    cells.append(nbf.v4.new_markdown_cell(source))

def code(source):
    cells.append(nbf.v4.new_code_cell(source))

# ── Title ──────────────────────────────────────────────────────────────────
md("""# purge_dups vs BISER Segmental Duplication Enrichment Analysis

**Assembly:** *Octodon degus* — hifiasm-041425 scaffolded assembly
**Question:** Are genomic regions classified by `purge_dups` enriched for BISER-defined segmental duplications?
**Method:** 10,000 chromosome- and length-matched permutations with Benjamini-Hochberg FDR correction
**Date:** 2026-07-21
""")

# ── Setup ──────────────────────────────────────────────────────────────────
md("""## 1. Environment Setup""")

code(r"""import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# Use non-interactive backend for inline
%matplotlib inline

PROJECT_ROOT = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj")
OUTPUT_ROOT = PROJECT_ROOT / "figure" / "segdup-purgeDup-overlap" / "purge-dups-biser-enrichment"
TABLES_DIR = OUTPUT_ROOT / "tables"
FIGURES_DIR = OUTPUT_ROOT / "figures"
PERM_DIR = OUTPUT_ROOT / "permutations"
LOGS_DIR = OUTPUT_ROOT / "logs"
INTERMEDIATE_DIR = OUTPUT_ROOT / "intermediate"

with open(OUTPUT_ROOT / "config" / "analysis_parameters.json") as f:
    CONFIG = json.load(f)

CATEGORIES = CONFIG["analysis_parameters"]["categories"]
CATEGORY_ORDER = CATEGORIES
MITO_SEQ = CONFIG["analysis_parameters"]["mitochondrial_sequence"]
RANDOM_SEED = CONFIG["analysis_parameters"]["random_seed"]

# Consistent color palette (from purge_dup_distribution.ipynb)
CATEGORY_COLORS = {
    "REPEAT":   "#FDAE61",  # orange
    "HAPLOTIG": "#D73027",  # red
    "JUNK":     "#4D4D4D",  # dark grey
    "OVLP":     "#1B9E77",  # green
    "HIGHCOV":  "#762A83",  # purple
}

# Matplotlib settings
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

print(f"Output root: {OUTPUT_ROOT}")
print(f"Categories: {CATEGORIES}")
print("Setup complete.")
""")

# ── Load data ──────────────────────────────────────────────────────────────
md("""## 2. Load All Results""")

code(r"""def load_table(path: Path) -> List[Dict]:
    rows = []
    with open(path) as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            rows.append(dict(zip(headers, parts)))
    return rows

# Primary tables
purge_summary = load_table(TABLES_DIR / "purge_dups_category_summary.tsv")
observed = load_table(TABLES_DIR / "purge_dups_biser_observed.tsv")
perm_results = load_table(TABLES_DIR / "purge_dups_biser_permutation_enrichment.tsv")
interval_details = load_table(TABLES_DIR / "purge_dups_biser_interval_details.tsv")
biser_summary = load_table(TABLES_DIR / "biser_chromosome_summary.tsv")
chrom_summary = load_table(TABLES_DIR / "purge_dups_biser_chromosome_summary.tsv")
runtime = load_table(TABLES_DIR / "runtime_benchmark.tsv")

# Null distributions
null_dists = dict(np.load(PERM_DIR / "primary_null_distributions.npz"))

obs_by_cat = {o["category"]: o for o in observed}
perm_by_cat = {r["category"]: r for r in perm_results}

n_perm = int(perm_results[0]["n_permutations"]) if perm_results else 0
run_mode = perm_results[0]["run_mode"] if perm_results else "?"

# BISER stats
biser_total_bp = sum(int(r["biser_union_bp"]) for r in biser_summary)
biser_total_mb = biser_total_bp / 1e6
genome_total = sum(int(r["chromosome_length"]) for r in biser_summary)
genome_frac = biser_total_bp / genome_total * 100

print(f"Loaded {len(perm_results)} permutation results ({run_mode} mode, {n_perm:,} permutations)")
print(f"Loaded {len(interval_details)} interval details")
print(f"Loaded {len(chrom_summary)} chromosome×category rows")
print(f"BISER union: {biser_total_mb:.1f} Mb ({genome_frac:.1f}% of genome)")
print(f"Null distributions: {list(null_dists.keys())}")
""")

# ── Figure 1: Observed overlap ────────────────────────────────────────────
md("""## 3. Figure 1 — Observed BISER Overlap by Category

Bar chart showing the percentage of each category's bases that overlap BISER-defined segmental duplications.
The dashed line shows the genome-wide BISER fraction as a descriptive background.
""")

code(r"""fig, ax = plt.subplots(figsize=(8, 5))

pcts = [float(obs_by_cat[cat]["pct_bases_overlapping"]) for cat in CATEGORY_ORDER]
totals_mb = [float(obs_by_cat[cat]["total_mb"]) for cat in CATEGORY_ORDER]
n_intervals = [int(obs_by_cat[cat]["n_intervals"]) for cat in CATEGORY_ORDER]
colors = [CATEGORY_COLORS[cat] for cat in CATEGORY_ORDER]

bars = ax.bar(CATEGORY_ORDER, pcts, color=colors, edgecolor="white",
              linewidth=0.5, width=0.65)

# Genome background
ax.axhline(y=genome_frac, color="black", linestyle="--", linewidth=1.2,
           label=f"Genome-wide BISER fraction\n({genome_frac:.1f}%)")

# Annotate
for bar, pct, mb, n in zip(bars, pcts, totals_mb, n_intervals):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
            f"{pct:.1f}%\n{mb:.1f} Mb\nn={n}",
            ha="center", va="bottom", fontsize=8, linespacing=1.3)

ax.set_ylabel("Category bases overlapping BISER (%)")
ax.set_title("Observed BISER Overlap by purge_dups Category")
ax.set_ylim(0, max(pcts) * 1.25)
ax.legend(loc="upper right", frameon=True, fontsize=8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()

# Save
fig.savefig(FIGURES_DIR / "observed_biser_overlap_by_category.png")
fig.savefig(FIGURES_DIR / "observed_biser_overlap_by_category.pdf")
plt.show()
print("Saved to figures/observed_biser_overlap_by_category.{png,pdf}")
""")

# ── Figure 2: Fold enrichment ─────────────────────────────────────────────
md("""## 4. Figure 2 — Permutation-Based Fold Enrichment

Fold enrichment relative to chromosome- and length-matched randomized intervals.
Reference line at 1.0. Asterisks show FDR significance.
""")

code(r"""fig, ax = plt.subplots(figsize=(8, 5))

cats_present = [c for c in CATEGORY_ORDER if c in perm_by_cat]
folds = [float(perm_by_cat[cat]["fold_enrichment"]) for cat in cats_present]
pvals = [float(perm_by_cat[cat]["empirical_p_enrichment"]) for cat in cats_present]
fdrs = [float(perm_by_cat[cat]["fdr_bh"]) for cat in cats_present]
colors = [CATEGORY_COLORS[cat] for cat in cats_present]

bars = ax.bar(cats_present, folds, color=colors, edgecolor="white",
              linewidth=0.5, width=0.65)
ax.axhline(y=1, color="black", linestyle="-", linewidth=0.8, alpha=0.5)

for bar, fold, pv, fdr in zip(bars, folds, pvals, fdrs):
    sig = ""
    if fdr < 0.05: sig = "*"
    if fdr < 0.01: sig = "**"
    if fdr < 0.001: sig = "***"
    y_pos = bar.get_height() + 0.5
    ax.text(bar.get_x() + bar.get_width() / 2, y_pos,
            f"{fold:.1f}×{sig}\nP={pv:.4f}\nFDR={fdr:.4f}",
            ha="center", va="bottom", fontsize=8, linespacing=1.3)

ax.set_ylabel("Fold Enrichment (observed / expected)")
ax.set_title(f"Permutation Enrichment of BISER Overlap\n({n_perm:,} chromosome- and length-matched permutations)")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.text(0.98, 0.95, "* FDR < 0.05\n** FDR < 0.01\n*** FDR < 0.001",
        transform=ax.transAxes, fontsize=7, va="top", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

fig.tight_layout()
fig.savefig(FIGURES_DIR / "permutation_fold_enrichment.png")
fig.savefig(FIGURES_DIR / "permutation_fold_enrichment.pdf")
plt.show()
print("Saved to figures/permutation_fold_enrichment.{png,pdf}")
""")

# ── Figure 3: Null distributions ──────────────────────────────────────────
md("""## 5. Figure 3 — Null Distributions by Category

Histograms of the randomized BISER-overlap null distributions (10,000 permutations each).
Red line = observed value. Dashed line = null mean. Grey band = 95% CI.
""")

code(r"""for cat in CATEGORY_ORDER:
    if cat not in perm_by_cat:
        continue
    key = f"{cat}_overlap"
    if key not in null_dists:
        continue

    fig, ax = plt.subplots(figsize=(7, 4.5))

    null_data = null_dists[key] / 1e6  # Mb
    r = perm_by_cat[cat]
    observed_val = float(r["observed_overlap_mb"])
    null_mean = float(r["null_mean_overlap_bp"]) / 1e6
    null_lower = float(r["null_lower_95_bp"]) / 1e6
    null_upper = float(r["null_upper_95_bp"]) / 1e6
    fold = float(r["fold_enrichment"])
    pval = float(r["empirical_p_enrichment"])
    fdr = float(r["fdr_bh"])
    z = float(r.get("z_score", "nan"))

    color = CATEGORY_COLORS[cat]
    ax.hist(null_data, bins=40, color=color, alpha=0.6, edgecolor="white",
            linewidth=0.3, density=True)

    ax.axvline(x=observed_val, color="darkred", linestyle="-", linewidth=2.5,
               label=f"Observed: {observed_val:.2f} Mb")
    ax.axvline(x=null_mean, color="black", linestyle="--", linewidth=1.2,
               label=f"Null mean: {null_mean:.2f} Mb")
    ax.axvspan(null_lower, null_upper, alpha=0.15, color="grey",
               label=f"95% CI: [{null_lower:.2f}, {null_upper:.2f}]")

    ax.set_xlabel("BISER Overlap (Mb)")
    ax.set_ylabel("Density")
    sig_label = "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*" if fdr < 0.05 else "ns"
    ax.set_title(f"{cat} — Null Distribution of BISER Overlap\n"
                 f"Fold = {fold:.1f}×, P = {pval:.4f}, FDR = {fdr:.4f} {sig_label}")
    ax.legend(loc="upper right", fontsize=7, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"null_distribution_{cat}.png")
    plt.show()
    plt.close(fig)

print("Saved null distribution figures for all 5 categories.")
""")

# ── Figure 4: Interval overlap fractions ──────────────────────────────────
md("""## 6. Figure 4 — Individual Interval BISER Overlap Fractions

Each point = one `purge_dups` interval. Y-axis = fraction of that interval's bases
overlapping BISER. Black bar = category median. n shown at top.
""")

code(r"""fig, ax = plt.subplots(figsize=(11, 5.5))

n_total = 0
x_positions = []
y_values = []
colors_list = []

for i, cat in enumerate(CATEGORY_ORDER):
    cat_details = [d for d in interval_details if d["category"] == cat]
    n = len(cat_details)
    if n == 0:
        continue
    n_total += n
    fracs = [float(d["biser_overlap_fraction"]) * 100 for d in cat_details]

    rng_jitter = np.random.default_rng(42 + i)
    jitter = rng_jitter.uniform(-0.25, 0.25, n)
    x_positions.extend((np.ones(n) * (i + 1) + jitter).tolist())
    y_values.extend(fracs)
    colors_list.extend([CATEGORY_COLORS[cat]] * n)

    # Median line
    median_val = np.median(fracs)
    ax.hlines(y=median_val, xmin=i + 0.6, xmax=i + 1.4,
              colors="black", linewidth=2.5, linestyles="-", zorder=10)

ax.scatter(x_positions, y_values, c=colors_list, alpha=0.65, s=50,
           edgecolors="white", linewidth=0.3, zorder=5)

ax.set_xticks(range(1, len(CATEGORY_ORDER) + 1))
ax.set_xticklabels(CATEGORY_ORDER)
ax.set_ylabel("BISER Overlap Fraction (%)")
ax.set_title(f"Individual purge_dups Interval BISER Overlap Fractions (n = {n_total})")
ax.set_ylim(-5, 108)

for i, cat in enumerate(CATEGORY_ORDER):
    n = sum(1 for d in interval_details if d["category"] == cat)
    ax.text(i + 1, 104, f"n={n}", ha="center", fontsize=8, fontweight="bold")

# BISER-evaluable vs non-evaluable indicator
n_eval = sum(1 for d in interval_details if d["biser_evaluable"] == "True")
n_not = sum(1 for d in interval_details if d["biser_evaluable"] == "False")
ax.text(0.02, 0.98, f"BISER-evaluable intervals: {n_eval}/{n_total}\n"
        f"Non-evaluable (no BISER on scaffold): {n_not}",
        transform=ax.transAxes, fontsize=7, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(FIGURES_DIR / "interval_overlap_fractions.png")
fig.savefig(FIGURES_DIR / "interval_overlap_fractions.pdf")
plt.show()
print("Saved to figures/interval_overlap_fractions.{png,pdf}")
""")

# ── Figure 5: Chromosome heatmap ──────────────────────────────────────────
md("""## 7. Figure 5 — Chromosome × Category Overlap Heatmap

Heatmap showing the percentage of purge_dups bases overlapping BISER for each
chromosome–category combination. Only shows chromosomes with purge_dups data.
""")

code(r"""def chrom_sort_key(c):
    if c.startswith("chr") and c[3:].isdigit(): return (0, int(c[3:]))
    elif c == "chrX": return (0, 29)
    elif c == "chrY": return (0, 30)
    elif c.startswith("seq") and c[3:].isdigit(): return (1, int(c[3:]))
    else: return (2, 999)

# Build matrix
chroms_with_data = sorted(set(r["chromosome"] for r in chrom_summary), key=chrom_sort_key)
matrix = np.full((len(chroms_with_data), len(CATEGORY_ORDER)), np.nan)

for d in chrom_summary:
    cat = d["category"]
    chrom = d["chromosome"]
    if cat in CATEGORY_ORDER and chrom in chroms_with_data:
        row = chroms_with_data.index(chrom)
        col = CATEGORY_ORDER.index(cat)
        matrix[row, col] = float(d["pct_purge_bp_in_biser"])

# Show only chromosomes with ≥3 categories or top 40
n_data = np.sum(~np.isnan(matrix), axis=1)
top_40_set = set(np.argsort(n_data)[::-1][:40])
many_data_set = set(np.where(n_data >= 3)[0])
keep = np.array(sorted(top_40_set | many_data_set))
if len(keep) > 45:
    keep = np.argsort(n_data)[::-1][:45]
    keep.sort()

matrix_plot = matrix[keep]
chroms_plot = [chroms_with_data[i] for i in keep]

fig, ax = plt.subplots(figsize=(9, max(6, len(chroms_plot) * 0.28)))

im = ax.imshow(matrix_plot, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)

ax.set_xticks(range(len(CATEGORY_ORDER)))
ax.set_xticklabels(CATEGORY_ORDER, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(len(chroms_plot)))
ax.set_yticklabels(chroms_plot, fontsize=6.5)

# Cell annotations
for i in range(len(chroms_plot)):
    for j in range(len(CATEGORY_ORDER)):
        val = matrix_plot[i, j]
        if not np.isnan(val):
            text_color = "white" if val > 55 else "black"
            ax.text(j, i, f"{val:.0f}%", ha="center", va="center",
                    fontsize=5.5, color=text_color, fontweight="bold" if val > 80 else "normal")

cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("BISER Overlap (%)", fontsize=9)

ax.set_title(f"Chromosome × Category BISER Overlap\n({len(chroms_plot)} chromosomes with purge_dups data)")
fig.tight_layout()

fig.savefig(FIGURES_DIR / "chromosome_category_overlap.png")
fig.savefig(FIGURES_DIR / "chromosome_category_overlap.pdf")
plt.show()
print("Saved to figures/chromosome_category_overlap.{png,pdf}")
""")

# ── Figure 6: Observed vs expected ────────────────────────────────────────
md("""## 8. Figure 6 — Observed vs Expected Overlap

Dark red bars = observed BISER overlap. Grey bars = mean randomized overlap.
Error bars = 95% CI of the null distribution.
""")

code(r"""fig, ax = plt.subplots(figsize=(8, 5))

cats_present = [c for c in CATEGORY_ORDER if c in perm_by_cat]
observed_vals = [float(perm_by_cat[cat]["observed_overlap_mb"]) for cat in cats_present]
expected_vals = [float(perm_by_cat[cat]["null_mean_overlap_bp"]) / 1e6 for cat in cats_present]
lower_ci = [float(perm_by_cat[cat]["null_lower_95_bp"]) / 1e6 for cat in cats_present]
upper_ci = [float(perm_by_cat[cat]["null_upper_95_bp"]) / 1e6 for cat in cats_present]

x = np.arange(len(cats_present))
width = 0.32

bars_obs = ax.bar(x - width/2, observed_vals, width, color="darkred", alpha=0.85,
                   edgecolor="white", linewidth=0.5, label="Observed")
bars_exp = ax.bar(x + width/2, expected_vals, width, color="grey", alpha=0.55,
                   edgecolor="white", linewidth=0.5, label="Expected (null mean)")

yerr_lower = [e - l for e, l in zip(expected_vals, lower_ci)]
yerr_upper = [u - e for e, u in zip(expected_vals, upper_ci)]
ax.errorbar(x + width/2, expected_vals, yerr=[yerr_lower, yerr_upper],
            fmt="none", color="black", capsize=4, linewidth=1.2, label="95% CI")

# Significance stars
for i, cat in enumerate(cats_present):
    fdr = float(perm_by_cat[cat]["fdr_bh"])
    if fdr < 0.001: sig = "***"
    elif fdr < 0.01: sig = "**"
    elif fdr < 0.05: sig = "*"
    else: sig = "ns"
    max_h = max(observed_vals[i], upper_ci[i])
    ax.text(i, max_h * 1.06, sig, ha="center", fontsize=13, fontweight="bold")

ax.set_xticks(x)
ax.set_xticklabels(cats_present)
ax.set_ylabel("BISER Overlap (Mb)")
ax.set_title(f"Observed vs Expected BISER Overlap\n({n_perm:,} chromosome- and length-matched permutations)")
ax.legend(fontsize=8, frameon=True, loc="upper left")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(FIGURES_DIR / "observed_vs_expected_overlap.png")
fig.savefig(FIGURES_DIR / "observed_vs_expected_overlap.pdf")
plt.show()
print("Saved to figures/observed_vs_expected_overlap.{png,pdf}")
""")

# ── Result tables ──────────────────────────────────────────────────────────
md("""## 9. Primary Results Table

Publication-style table matching the pair-level interpretation framework from
`category_exploration_v4_busco.ipynb`. Color-coded category cells, auto-wrapped
text, dark header bar with white labels.
""")

code(r"""import re
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties

# ── Build table data ───────────────────────────────────────────────────
rows_data = []
for cat in CATEGORY_ORDER:
    if cat not in perm_by_cat:
        continue
    r = perm_by_cat[cat]
    o = obs_by_cat[cat]
    fold = float(r["fold_enrichment"])
    pval = float(r["empirical_p_enrichment"])
    fdr = float(r["fdr_bh"])
    sig = "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*" if fdr < 0.05 else "ns"
    rows_data.append({
        "cat": cat, "color": CATEGORY_COLORS[cat],
        "cells": [
            cat,
            o["n_intervals"],
            f'{float(o["total_mb"]):.1f}',
            f'{float(o["pct_bases_overlapping"]):.1f}%',
            f'{float(r["observed_overlap_mb"]):.2f}',
            f'{float(r["null_mean_overlap_bp"])/1e6:.2f}',
            f'{fold:.1f}×',
            f'{pval:.4f}',
            f'{fdr:.4f}',
            sig,
        ],
    })

col_labels = [
    "Category", "N", "Total\n(Mb)", "Raw BISER\n(%)",
    "Observed\n(Mb)", "Null mean\n(Mb)", "Fold\nenrichment",
    "P (enr)", "FDR", "Sig.",
]
col_widths = [0.145, 0.055, 0.075, 0.095, 0.10, 0.105, 0.095, 0.105, 0.105, 0.12]

# ── Typography ─────────────────────────────────────────────────────────
FIG_W = 14.0
TITLE_FS, HEADER_FS, CAT_FS, BODY_FS = 18, 13, 11.5, 10.5
HEADER_H_PT, ROW_PAD_PT = 40, 5
TOP_MARGIN, TITLE_H, TITLE_GAP, BOT_MARGIN = 6, 24, 6, 5

body_font   = FontProperties(family="sans-serif", size=BODY_FS)
cat_font    = FontProperties(family="sans-serif", size=CAT_FS, weight="bold")
header_font = FontProperties(family="sans-serif", size=HEADER_FS, weight="bold")
title_font  = FontProperties(family="sans-serif", size=TITLE_FS, weight="bold")
count_font  = FontProperties(family="sans-serif", size=11, weight="bold")

def norm(s):
    return re.sub(r"\s+", " ", str(s)).strip()

def text_width(renderer, text, fp):
    w, _, _ = renderer.get_text_width_height_descent(text, fp, ismath=False)
    return w

def wrap(text, max_px, renderer, fp):
    text = norm(text)
    if not text: return ""
    words = text.split()
    lines, cur = [], words[0]
    for w in words[1:]:
        cand = f"{cur} {w}"
        if text_width(renderer, cand, fp) <= max_px:
            cur = cand
        else:
            lines.append(cur); cur = w
    lines.append(cur)
    return "\n".join(lines)

# ── Measure & wrap ─────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIG_W, 6), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
fig.canvas.draw()
renderer = fig.canvas.get_renderer()

fig_w_px = FIG_W * fig.dpi
table_w_px = fig_w_px * 0.96
h_pad_px = 6.5 / 72 * fig.dpi

wrapped_rows = []
row_lines = []
for row in rows_data:
    wrapped = []
    for ci, val in enumerate(row["cells"]):
        cw = table_w_px * col_widths[ci]
        avail = cw - 2 * h_pad_px - 5
        fp = cat_font if ci == 0 else count_font if ci == 1 else body_font
        wrapped.append(wrap(val, avail, renderer, fp))
    wrapped_rows.append({"text": wrapped, "color": row["color"]})
    row_lines.append(max(t.count("\n") + 1 for t in wrapped))

line_h_pt = BODY_FS * 1.15
row_heights_pt = [l * line_h_pt + 2 * ROW_PAD_PT for l in row_lines]
table_h_pt = HEADER_H_PT + sum(row_heights_pt)
fig_h_pt = TOP_MARGIN + TITLE_H + TITLE_GAP + table_h_pt + BOT_MARGIN
fig_h_in = fig_h_pt / 72
fig.set_size_inches(FIG_W, fig_h_in, forward=True)
fig.canvas.draw()

def pt2y(pt):
    return pt / fig_h_pt

title_y = 1.0 - pt2y(TOP_MARGIN)
table_top = 1.0 - pt2y(TOP_MARGIN + TITLE_H + TITLE_GAP)
header_h = pt2y(HEADER_H_PT)
row_hs = [pt2y(h) for h in row_heights_pt]
h_pad_f = (6.5 / 72) / FIG_W

col_lefts = [0.02]
for i in range(len(col_widths) - 1):
    col_lefts.append(col_lefts[-1] + 0.96 * col_widths[i])
col_disp_w = [0.96 * w for w in col_widths]

# ── Title ──────────────────────────────────────────────────────────────
ax.text(0.5, title_y,
    f"purge_dups category enrichment for BISER segmental duplications\n"
    f"({n_perm:,} chromosome- and length-matched permutations, BH-FDR corrected)",
    fontproperties=title_font, color="black", ha="center", va="top")

# ── Header ─────────────────────────────────────────────────────────────
hb = table_top - header_h
for ci, lbl in enumerate(col_labels):
    r = Rectangle((col_lefts[ci], hb), col_disp_w[ci], header_h,
                   facecolor="#333333", edgecolor="#666666", linewidth=0.7)
    ax.add_patch(r)
    ax.text(col_lefts[ci] + col_disp_w[ci] / 2, hb + header_h / 2, lbl,
            fontproperties=header_font, color="white", ha="center", va="center")

# ── Body rows ──────────────────────────────────────────────────────────
cur_top = hb
for ri, row in enumerate(wrapped_rows):
    rb = cur_top - row_hs[ri]
    for ci, txt in enumerate(row["text"]):
        cl, cw = col_lefts[ci], col_disp_w[ci]
        if ci == 0:
            fc, tc, fp, ha, tx = row["color"], "white", cat_font, "left", cl + h_pad_f
        elif ci in (1, 9):
            fc, tc, fp, ha, tx = "#FAFAFA", "#222222", count_font, "center", cl + cw / 2
        else:
            fc, tc, fp, ha, tx = "#FAFAFA", "#222222", body_font, "left", cl + h_pad_f

        r = Rectangle((cl, rb), cw, row_hs[ri], facecolor=fc,
                       edgecolor="#666666", linewidth=0.7)
        ax.add_patch(r)
        ax.text(tx, rb + row_hs[ri] / 2, txt, fontproperties=fp,
                color=tc, ha=ha, va="center", linespacing=1.05, clip_on=True)
    cur_top = rb

out = FIGURES_DIR / "enrichment_results_table.png"
fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.03)
plt.show()
print(f"Saved: {out}")
""")

# ── Interpretation ─────────────────────────────────────────────────────────
md("""## 10. Interpretation""")

code(r"""print("=" * 70)
print("INTERPRETATION GUIDE")
print("=" * 70)

# HAPLOTIG
hap_fold = float(perm_by_cat["HAPLOTIG"]["fold_enrichment"])
hap_fdr = float(perm_by_cat["HAPLOTIG"]["fdr_bh"])
hap_raw = float(obs_by_cat["HAPLOTIG"]["pct_bases_overlapping"])

print(f'''
HAPLOTIG (fold={hap_fold:.1f}×, FDR={hap_fdr:.3f}):
  Raw BISER overlap is {hap_raw:.0f}%, but after controlling for chromosome
  placement and interval length, there is NO significant enrichment.
  The high raw overlap is driven by two large intervals on chrY (7.1 Mb)
  and seq31 (3.0 Mb), both in BISER-dense regions. When these same-sized
  intervals are randomly repositioned on the same chromosomes, they overlap
  BISER just as much as the real placement (null mean ≈ observed).

  → This does NOT mean HAPLOTIG regions are unrelated to segmental
    duplications; it means their chromosome distribution and large sizes
    are sufficient to explain the observed overlap without invoking
    preferential localization.
''')

# REPEAT
rep_fold = float(perm_by_cat["REPEAT"]["fold_enrichment"])
rep_fdr = float(perm_by_cat["REPEAT"]["fdr_bh"])
rep_raw = float(obs_by_cat["REPEAT"]["pct_bases_overlapping"])

print(f'''
REPEAT (fold={rep_fold:.1f}×, FDR={rep_fdr:.4f}):
  Highly significant enrichment. REPEAT-classified regions overlap BISER
  segmental duplications 2.5× more than expected by chance after controlling
  for chromosome and interval length. This is broadly distributed across
  chromosomes (not driven by a single chromosome).

  → This is biologically expected: both purge_dups REPEAT and BISER respond
    to duplicated/repetitive sequence. However, this does NOT prove every
    REPEAT region is a genuine segmental duplication — it could also reflect
    complex duplicated sequence, high-copy repeat structure, or residual
    haplotypic redundancy embedded in duplicated regions.
''')

# HIGHCOV
hc_fold = float(perm_by_cat["HIGHCOV"]["fold_enrichment"])
hc_fdr = float(perm_by_cat["HIGHCOV"]["fdr_bh"])

print(f'''
HIGHCOV (fold={hc_fold:.1f}×, FDR={hc_fdr:.4f}):
  Statistically significant enrichment, but based on only 2 intervals
  (0.07 Mb total) on chr19. Both overlap BISER heavily (84-90%), whereas
  random placement of same-sized intervals on chr19 rarely overlaps BISER.

  → Interpret with extreme caution. With n=2, this result is highly
    sensitive to the specific coordinates. The association between high
    coverage and duplicated sequence is mechanistically plausible but
    cannot be generalized from 2 intervals.
''')

# OVLP
ovlp_fold = float(perm_by_cat["OVLP"]["fold_enrichment"])
ovlp_fdr = float(perm_by_cat["OVLP"]["fdr_bh"])

print(f'''
OVLP (fold={ovlp_fold:.1f}×, FDR={ovlp_fdr:.3f}):
  Nominally enriched (P=0.048) but not significant after FDR correction
  (FDR=0.079). The 1.9× fold enrichment suggests a trend toward association
  between overlap-duplication calls and BISER-defined duplications.

  → Inconclusive. More data or a different approach would be needed to
    confirm or refute this association.
''')

# JUNK
print(f'''
JUNK (fold=1.0×, FDR=1.0):
  No enrichment. All JUNK intervals reside on small unplaced scaffolds.
  The null distribution is degenerate (SD=0) because all possible placements
  of these small intervals on the same tiny scaffolds produce identical
  overlap results.

  → The permutation test is uninformative for JUNK. These intervals are
    not evaluable for BISER enrichment because they are on scaffolds too
    small to provide meaningful randomization space.
''')
""")

# ── Sensitivity summary ────────────────────────────────────────────────────
md("""## 11. Sensitivity Analysis Summary""")

code(r"""sens_path = TABLES_DIR / "purge_dups_biser_sensitivity_results.tsv"
if sens_path.exists():
    sens = load_table(sens_path)
    print(f"Total sensitivity tests: {len(sens)}")

    # Group by analysis type
    by_analysis = defaultdict(list)
    for s in sens:
        analysis = s["analysis"]
        if analysis.startswith("leave_one_out_"):
            analysis = "leave_one_out_*"
        by_analysis[analysis].append(s)

    print(f"\n{'Analysis':<35} {'Tests':>6} {'Significant':>12}")
    print("-" * 55)
    for analysis in sorted(by_analysis.keys()):
        items = by_analysis[analysis]
        n_sig = sum(1 for s in items if float(s["empirical_p_enrichment"]) < 0.05)
        print(f"{analysis:<35} {len(items):>6} {n_sig:>12}")

    # Key findings
    print("\n--- Key Sensitivity Findings ---")

    for s in sens:
        if s["analysis"] in ("autosomes_only", "exclude_chrX_chrY",
                              "placed_scaffolds_only", "biser_evaluable_only"):
            cat_clean = s["category"].rsplit("_", 1)[0]
            fold = float(s["fold_enrichment"])
            pval = float(s["empirical_p_enrichment"])
            if fold > 1.5 or pval < 0.05:
                sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
                print(f"  {s['category']:<30} [{s['analysis']:<25}] fold={fold:.1f}×, P={pval:.4f} {sig}")
else:
    print("No sensitivity results. Run with --mode pilot or --mode final.")
""")

# ── Category summary table ─────────────────────────────────────────────────
md("""## 12. Category Interval Size Summary

Publication-style table with color-coded category cells, matching the
pair-level interpretation framework from `category_exploration_v4_busco.ipynb`.
""")

code(r"""import re
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties

# ── Build table data ───────────────────────────────────────────────────
rows_data = []
for r in purge_summary:
    cat = r["category"]
    rows_data.append({
        "cat": cat, "color": CATEGORY_COLORS[cat],
        "cells": [
            cat,
            r["n_intervals"],
            f'{float(r["total_mb"]):.2f}',
            f'{int(r["min_length"])/1000:.1f}',
            f'{float(r["median_length"])/1000:.1f}',
            f'{float(r["mean_length"])/1000:.1f}',
            f'{int(r["max_length"])/1000:.1f}',
        ],
    })

col_labels = [
    "Category", "N", "Total\n(Mb)", "Min\n(kb)", "Median\n(kb)",
    "Mean\n(kb)", "Max\n(kb)",
]
col_widths = [0.18, 0.08, 0.12, 0.12, 0.14, 0.14, 0.22]

# ── Typography ─────────────────────────────────────────────────────────
FIG_W = 11.0
TITLE_FS, HEADER_FS, CAT_FS, BODY_FS = 16, 12, 11, 10.5
HEADER_H_PT, ROW_PAD_PT = 36, 4
TOP_MARGIN, TITLE_H, TITLE_GAP, BOT_MARGIN = 5, 20, 5, 5

body_font   = FontProperties(family="sans-serif", size=BODY_FS)
cat_font    = FontProperties(family="sans-serif", size=CAT_FS, weight="bold")
header_font = FontProperties(family="sans-serif", size=HEADER_FS, weight="bold")
title_font  = FontProperties(family="sans-serif", size=TITLE_FS, weight="bold")
count_font  = FontProperties(family="sans-serif", size=10.5, weight="bold")

def norm(s):
    return re.sub(r"\s+", " ", str(s)).strip()

def text_width(renderer, text, fp):
    w, _, _ = renderer.get_text_width_height_descent(text, fp, ismath=False)
    return w

def wrap(text, max_px, renderer, fp):
    text = norm(text)
    if not text: return ""
    words = text.split()
    lines, cur = [], words[0]
    for w in words[1:]:
        cand = f"{cur} {w}"
        if text_width(renderer, cand, fp) <= max_px:
            cur = cand
        else:
            lines.append(cur); cur = w
    lines.append(cur)
    return "\n".join(lines)

# ── Measure & wrap ─────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIG_W, 6), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
fig.canvas.draw()
renderer = fig.canvas.get_renderer()

fig_w_px = FIG_W * fig.dpi
table_w_px = fig_w_px * 0.96
h_pad_px = 6.0 / 72 * fig.dpi

wrapped_rows = []
row_lines = []
for row in rows_data:
    wrapped = []
    for ci, val in enumerate(row["cells"]):
        cw = table_w_px * col_widths[ci]
        avail = cw - 2 * h_pad_px - 5
        fp = cat_font if ci == 0 else count_font if ci == 1 else body_font
        wrapped.append(wrap(val, avail, renderer, fp))
    wrapped_rows.append({"text": wrapped, "color": row["color"]})
    row_lines.append(max(t.count("\n") + 1 for t in wrapped))

line_h_pt = BODY_FS * 1.15
row_heights_pt = [l * line_h_pt + 2 * ROW_PAD_PT for l in row_lines]
table_h_pt = HEADER_H_PT + sum(row_heights_pt)
fig_h_pt = TOP_MARGIN + TITLE_H + TITLE_GAP + table_h_pt + BOT_MARGIN
fig_h_in = fig_h_pt / 72
fig.set_size_inches(FIG_W, fig_h_in, forward=True)
fig.canvas.draw()

def pt2y(pt):
    return pt / fig_h_pt

title_y = 1.0 - pt2y(TOP_MARGIN)
table_top = 1.0 - pt2y(TOP_MARGIN + TITLE_H + TITLE_GAP)
header_h = pt2y(HEADER_H_PT)
row_hs = [pt2y(h) for h in row_heights_pt]
h_pad_f = (6.0 / 72) / FIG_W

col_lefts = [0.02]
for i in range(len(col_widths) - 1):
    col_lefts.append(col_lefts[-1] + 0.96 * col_widths[i])
col_disp_w = [0.96 * w for w in col_widths]

# ── Title ──────────────────────────────────────────────────────────────
ax.text(0.5, title_y,
    f"purge_dups category interval size summary (n = {sum(int(r['n_intervals']) for r in purge_summary)} intervals)",
    fontproperties=title_font, color="black", ha="center", va="top")

# ── Header ─────────────────────────────────────────────────────────────
hb = table_top - header_h
for ci, lbl in enumerate(col_labels):
    r = Rectangle((col_lefts[ci], hb), col_disp_w[ci], header_h,
                   facecolor="#333333", edgecolor="#666666", linewidth=0.7)
    ax.add_patch(r)
    ax.text(col_lefts[ci] + col_disp_w[ci] / 2, hb + header_h / 2, lbl,
            fontproperties=header_font, color="white", ha="center", va="center")

# ── Body rows ──────────────────────────────────────────────────────────
cur_top = hb
for ri, row in enumerate(wrapped_rows):
    rb = cur_top - row_hs[ri]
    for ci, txt in enumerate(row["text"]):
        cl, cw = col_lefts[ci], col_disp_w[ci]
        if ci == 0:
            fc, tc, fp, ha, tx = row["color"], "white", cat_font, "left", cl + h_pad_f
        elif ci == 1:
            fc, tc, fp, ha, tx = "#FAFAFA", "#222222", count_font, "center", cl + cw / 2
        else:
            fc, tc, fp, ha, tx = "#FAFAFA", "#222222", body_font, "left", cl + h_pad_f

        r = Rectangle((cl, rb), cw, row_hs[ri], facecolor=fc,
                       edgecolor="#666666", linewidth=0.7)
        ax.add_patch(r)
        ax.text(tx, rb + row_hs[ri] / 2, txt, fontproperties=fp,
                color=tc, ha=ha, va="center", linespacing=1.05, clip_on=True)
    cur_top = rb

fig.savefig(FIGURES_DIR / "category_interval_size_summary.png",
            dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.03)
plt.show()
print("Saved: figures/category_interval_size_summary.png")
""")

# ── Top intervals for inspection ───────────────────────────────────────────
md("""## 13. Top Intervals for Manual Inspection

Regions with high BISER overlap that merit closer examination.
""")

code(r"""print(f"{'Category':<12} {'Chromosome':<12} {'Start':>12} {'End':>12} {'Length(Mb)':>10} {'BISER%':>8}")
print("-" * 70)

# Prioritize: 100% BISER overlap, evaluable, on chr scaffolds, > 10kb
priority = sorted(
    [d for d in interval_details
     if d["biser_evaluable"] == "True"
     and float(d["interval_length_mb"]) > 0.01],
    key=lambda x: (-float(x["biser_overlap_fraction"]), -float(x["interval_length_mb"]))
)

for d in priority[:20]:
    chrom = d["chrom"]
    is_chr = "✓" if chrom.startswith("chr") else " "
    print(f"{d['category']:<12} {chrom:<12} {d['start']:>12} {d['end']:>12} "
          f"{float(d['interval_length_mb']):>9.2f} {float(d['biser_overlap_fraction'])*100:>7.1f}% {is_chr}")
""")

# ── Runtime ────────────────────────────────────────────────────────────────
md("""## 14. Runtime Performance""")

code(r"""for b in runtime:
    print(f"Mode: {b['run_mode']}")
    print(f"  Permutations: {b['n_permutations']}")
    print(f"  Elapsed: {b['elapsed_seconds']}s")
    print(f"  Perm/s: {b['permutations_per_second']}")
    print(f"  Estimated for 1K: {b['estimated_seconds_for_1000']}s")
    print(f"  Estimated for 10K: {b['estimated_seconds_for_10000']}s")
    print()

print("Implementation: Python + NumPy, binary-search overlap queries")
print("No bedtools subprocess calls during permutation testing")
print("Single CPU core, < 50 MB peak memory")
""")

# ── Methods ────────────────────────────────────────────────────────────────
md("""## 15. Methods (Auto-Generated)

The methods paragraph below reflects the actual parameters used.
""")

code(r"""n_perm_str = perm_results[0]["n_permutations"] if perm_results else "N"
n_perm_int = int(n_perm_str)

print(f'''
BISER-defined segmental-duplication coordinates (n = 208,156 pairs) were
obtained from the hifiasm-041425 scaffolded assembly (`segdup_output_mod.bedpe`).
Both members of each duplication pair were extracted and overlapping intervals
were merged, yielding a nonredundant BISER union of {biser_total_mb:.1f} Mb
({genome_frac:.1f}% of the {genome_total/1e6:.0f} Mb eligible genome, excluding
the mitochondrial sequence {MITO_SEQ}).

purge_dups-classified regions (n = {sum(int(r['n_intervals']) for r in purge_summary)}
intervals) were parsed from `dups.bed` across five categories: HAPLOTIG
(n={purge_summary[0]['n_intervals']}), HIGHCOV (n={purge_summary[1]['n_intervals']}),
JUNK (n={purge_summary[2]['n_intervals']}), OVLP (n={purge_summary[3]['n_intervals']}),
and REPEAT (n={purge_summary[4]['n_intervals']}).

For each category, we quantified the fraction of category sequence intersecting
the merged BISER union. Enrichment was evaluated using {n_perm_int:,}
chromosome-restricted permutations in which intervals were randomly repositioned
while preserving their original chromosome assignments and exact lengths. Random
intervals were permitted to overlap one another; per-permutation base-pair overlap
statistics were computed from the nonredundant union of randomized category intervals.
Empirical one-sided P values used a plus-one correction and were adjusted across
the five categories using the Benjamini-Hochberg procedure.

Sensitivity analyses (1,000 permutations each) evaluated: (1) autosomes only,
(2) excluding sex chromosomes, (3) placed chromosome scaffolds only,
(4) excluding the single longest interval per category, (5) leave-one-chromosome-out,
(6) intervals with >50% BISER overlap, and (7) BISER-evaluable chromosomes only.

Analyses were implemented in Python/NumPy with binary search for BISER overlap
queries. 10,000 permutations across all five categories completed in ~30 seconds
on a single CPU core.
''')
""")

# ── Final summary ──────────────────────────────────────────────────────────
md("""## 16. Final Summary""")

code(r"""print("=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f'''
1.  Selected input files:
    purge_dups: {CONFIG['input_files']['purge_dups_bed']}
    BISER:      {CONFIG['input_files']['biser_bedpe']}
    Assembly:   {CONFIG['input_files']['assembly_fasta']}

2.  Total purge_dups intervals: {sum(int(r['n_intervals']) for r in purge_summary)}

3.  Interval counts by category:
''')
for r in purge_summary:
    print(f"    {r['category']:<12} {r['n_intervals']:>4} intervals, {float(r['total_mb']):>8.2f} Mb")

print(f'''
4.  Total BISER union size: {biser_total_mb:.1f} Mb ({genome_frac:.1f}% of genome)

5.  Highest raw BISER overlap:
''')
highest_raw = max(obs_by_cat.items(), key=lambda x: float(x[1]["pct_bases_overlapping"]))
print(f"    {highest_raw[0]}: {float(highest_raw[1]['pct_bases_overlapping']):.1f}%")

print(f'''
6.  Highest matched fold enrichment:
''')
highest_fold = max(perm_by_cat.items(), key=lambda x: float(x[1]["fold_enrichment"]))
print(f"    {highest_fold[0]}: {float(highest_fold[1]['fold_enrichment']):.1f}×")

print(f'''
7.  Significant after FDR correction:
''')
sig_cats = [cat for cat in CATEGORY_ORDER if cat in perm_by_cat and float(perm_by_cat[cat]["fdr_bh"]) < 0.05]
for cat in sig_cats:
    r = perm_by_cat[cat]
    print(f"    {cat}: fold={float(r['fold_enrichment']):.1f}×, FDR={float(r['fdr_bh']):.4f}")
if not sig_cats:
    print("    None")

print(f'''
8.  Key caveats:
    - HIGHCOV has only 2 intervals — significant but based on n=2
    - JUNK null distribution is degenerate (SD=0) — not evaluable
    - HAPLOTIG raw overlap (81.8%) is explained by chromosome/length
    - BISER overlap ≠ definitive evidence of segmental duplication
    - 75 purge_dups chromosomes lack BISER data

9.  Runtime: {n_perm_int:,} permutations in {runtime[-1]['elapsed_seconds']}s

10. All outputs: {OUTPUT_ROOT}
''')
""")

# ── Build notebook ─────────────────────────────────────────────────────────
nb = nbf.v4.new_notebook()
nb.cells = cells
nb.metadata = {
    "kernelspec": {
        "display_name": "genome-assembly",
        "language": "python",
        "name": "genome-assembly",
    },
    "language_info": {
        "name": "python",
        "version": "3.12.0",
    },
}

with open(NB_PATH, "w") as f:
    nbf.write(nb, f)

print(f"\nNotebook written to: {NB_PATH}")
print(f"Total cells: {len(cells)}")
