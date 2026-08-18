#!/usr/bin/env python3
"""
03_plot_results.py — Generate publication-ready figures for purge_dups vs BISER enrichment.

Reads cached permutation results and summary tables. Does NOT recompute permutations.

Usage:
    python scripts/03_plot_results.py
"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj")
OUTPUT_ROOT = PROJECT_ROOT / "figure" / "segdup-purgeDup-overlap" / "purge-dups-biser-enrichment"
CONFIG_DIR = OUTPUT_ROOT / "config"
TABLES_DIR = OUTPUT_ROOT / "tables"
FIGURES_DIR = OUTPUT_ROOT / "figures"
PERM_DIR = OUTPUT_ROOT / "permutations"

with open(CONFIG_DIR / "analysis_parameters.json") as f:
    CONFIG = json.load(f)

CATEGORIES = CONFIG["analysis_parameters"]["categories"]
CATEGORY_ORDER = CATEGORIES

# ── Consistent colors ─────────────────────────────────────────────────────
CATEGORY_COLORS = {
    "REPEAT":   "#FDAE61",  # orange
    "HAPLOTIG": "#D73027",  # red
    "JUNK":     "#4D4D4D",  # dark grey
    "OVLP":     "#1B9E77",  # green
    "HIGHCOV":  "#762A83",  # purple
}

# ── Matplotlib settings ───────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
})


def load_observed_table() -> Dict[str, Dict]:
    """Load observed overlap table."""
    path = TABLES_DIR / "purge_dups_biser_observed.tsv"
    data: Dict[str, Dict] = {}
    with open(path) as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(headers, parts))
            data[row["category"]] = row
    return data


def load_permutation_table() -> List[Dict]:
    """Load permutation enrichment table."""
    path = TABLES_DIR / "purge_dups_biser_permutation_enrichment.tsv"
    results: List[Dict] = []
    with open(path) as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(headers, parts))
            results.append(row)
    return results


def load_null_distributions() -> Dict[str, np.ndarray]:
    """Load cached null distributions."""
    path = PERM_DIR / "primary_null_distributions.npz"
    if not path.exists():
        print(f"WARNING: Null distributions not found at {path}")
        return {}
    return dict(np.load(path))


def load_chromosome_table() -> List[Dict]:
    """Load chromosome-level summary."""
    path = TABLES_DIR / "purge_dups_biser_chromosome_summary.tsv"
    results: List[Dict] = []
    with open(path) as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(headers, parts))
            results.append(row)
    return results


def load_interval_details() -> List[Dict]:
    """Load interval-level details."""
    path = TABLES_DIR / "purge_dups_biser_interval_details.tsv"
    results: List[Dict] = []
    with open(path) as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(headers, parts))
            results.append(row)
    return results


# ── Figure 1: Observed BISER overlap ───────────────────────────────────────

def figure1_observed_overlap(obs: Dict[str, Dict]) -> str:
    """Bar chart: percentage of category bases overlapping BISER."""
    fig, ax = plt.subplots(figsize=(7, 5))

    categories = CATEGORY_ORDER
    pcts = [float(obs[cat]["pct_bases_overlapping"]) for cat in categories]
    totals_mb = [float(obs[cat]["total_mb"]) for cat in categories]
    n_intervals = [int(obs[cat]["n_intervals"]) for cat in categories]
    colors = [CATEGORY_COLORS[cat] for cat in categories]

    bars = ax.bar(categories, pcts, color=colors, edgecolor="white", linewidth=0.5, width=0.65)

    # Genome background line
    genome_frac = float(obs[categories[0]]["eligible_genome_biser_frac"]) * 100
    ax.axhline(y=genome_frac, color="black", linestyle="--", linewidth=1.2,
               label=f"Genome-wide BISER\n({genome_frac:.1f}%)")

    # Annotate bars
    for bar, pct, mb, n in zip(bars, pcts, totals_mb, n_intervals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
                f"{pct:.1f}%\n{mb:.1f} Mb\nn={n}",
                ha="center", va="bottom", fontsize=8, linespacing=1.3)

    ax.set_ylabel("Category bases overlapping BISER (%)")
    ax.set_title("Observed BISER Overlap by purge_dups Category")
    ax.set_ylim(0, max(pcts) * 1.25)
    ax.legend(loc="upper right", frameon=True, fontsize=8)

    # Remove top/right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    png_path = FIGURES_DIR / "observed_biser_overlap_by_category.png"
    pdf_path = FIGURES_DIR / "observed_biser_overlap_by_category.pdf"
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"  Figure 1: {png_path}")
    return str(png_path)


# ── Figure 2: Permutation fold enrichment ──────────────────────────────────

def figure2_fold_enrichment(perm_results: List[Dict]) -> str:
    """Bar chart: fold enrichment with significance annotations."""
    fig, ax = plt.subplots(figsize=(7, 5))

    perm_by_cat = {r["category"]: r for r in perm_results}
    categories = [c for c in CATEGORY_ORDER if c in perm_by_cat]

    folds = [float(perm_by_cat[cat]["fold_enrichment"]) for cat in categories]
    pvals = [float(perm_by_cat[cat]["empirical_p_enrichment"]) for cat in categories]
    fdrs = [float(perm_by_cat[cat]["fdr_bh"]) for cat in categories]
    colors = [CATEGORY_COLORS[cat] for cat in categories]

    bars = ax.bar(categories, folds, color=colors, edgecolor="white", linewidth=0.5, width=0.65)

    # Reference line at 1
    ax.axhline(y=1, color="black", linestyle="-", linewidth=0.8, alpha=0.5)

    # Annotate bars with p-value and significance
    for bar, fold, pv, fdr in zip(bars, folds, pvals, fdrs):
        sig = ""
        if fdr < 0.05:
            sig = "*"
        if fdr < 0.01:
            sig = "**"
        if fdr < 0.001:
            sig = "***"
        y_pos = bar.get_height() + 0.3
        ax.text(bar.get_x() + bar.get_width() / 2, y_pos,
                f"{fold:.1f}×{sig}\nP={pv:.4f}\nFDR={fdr:.4f}",
                ha="center", va="bottom", fontsize=8, linespacing=1.3)

    ax.set_ylabel("Fold Enrichment (observed / expected)")
    ax.set_title("Permutation-Based BISER Enrichment by purge_dups Category")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Add note about significance
    ax.text(0.98, 0.95, "* FDR < 0.05\n** FDR < 0.01\n*** FDR < 0.001",
            transform=ax.transAxes, fontsize=7, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    fig.tight_layout()

    png_path = FIGURES_DIR / "permutation_fold_enrichment.png"
    pdf_path = FIGURES_DIR / "permutation_fold_enrichment.pdf"
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"  Figure 2: {png_path}")
    return str(png_path)


# ── Figure 3: Null distributions ──────────────────────────────────────────

def figure3_null_distributions(
    perm_results: List[Dict],
    null_dists: Dict[str, np.ndarray],
) -> List[str]:
    """Histograms of null distributions per category with observed lines."""
    paths = []
    perm_by_cat = {r["category"]: r for r in perm_results}

    for cat in CATEGORY_ORDER:
        if cat not in perm_by_cat:
            continue
        key = f"{cat}_overlap"
        if key not in null_dists:
            continue

        fig, ax = plt.subplots(figsize=(6, 4.5))

        null_data = null_dists[key] / 1e6  # Convert to Mb
        r = perm_by_cat[cat]
        observed = float(r["observed_overlap_mb"])
        null_mean = float(r["null_mean_overlap_bp"]) / 1e6
        null_lower = float(r["null_lower_95_bp"]) / 1e6
        null_upper = float(r["null_upper_95_bp"]) / 1e6
        fold = float(r["fold_enrichment"])
        pval = float(r["empirical_p_enrichment"])
        fdr = float(r["fdr_bh"])

        color = CATEGORY_COLORS[cat]
        ax.hist(null_data, bins=30, color=color, alpha=0.6, edgecolor="white",
                linewidth=0.3, density=True)

        # Observed line
        ax.axvline(x=observed, color="darkred", linestyle="-", linewidth=2,
                   label=f"Observed: {observed:.2f} Mb")

        # Null mean
        ax.axvline(x=null_mean, color="black", linestyle="--", linewidth=1,
                   label=f"Null mean: {null_mean:.2f} Mb")

        # 95% CI
        ax.axvspan(null_lower, null_upper, alpha=0.15, color="grey",
                   label=f"95% CI: [{null_lower:.2f}, {null_upper:.2f}]")

        ax.set_xlabel("BISER Overlap (Mb)")
        ax.set_ylabel("Density")
        ax.set_title(f"{cat}: Null Distribution of BISER Overlap\n"
                     f"Fold={fold:.1f}×, P={pval:.4f}, FDR={fdr:.4f}")
        ax.legend(loc="upper right", fontsize=7, frameon=True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        fig.tight_layout()

        png_path = FIGURES_DIR / f"null_distribution_{cat}.png"
        fig.savefig(png_path)
        plt.close(fig)
        paths.append(str(png_path))
        print(f"  Figure 3 ({cat}): {png_path}")

    return paths


# ── Figure 4: Interval overlap fractions ──────────────────────────────────

def figure4_interval_overlap_fractions(interval_details: List[Dict]) -> str:
    """Scatter/strip plot of individual interval BISER overlap fractions."""
    fig, ax = plt.subplots(figsize=(10, 5))

    # Group by category
    by_cat: Dict[str, List[float]] = {c: [] for c in CATEGORY_ORDER}
    by_cat_chrom: Dict[str, List[str]] = {c: [] for c in CATEGORY_ORDER}
    for d in interval_details:
        cat = d["category"]
        if cat in by_cat:
            by_cat[cat].append(float(d["biser_overlap_fraction"]) * 100)

    # Plot as individual points with jitter
    x_positions = []
    y_values = []
    colors_list = []
    for i, cat in enumerate(CATEGORY_ORDER):
        n = len(by_cat[cat])
        if n == 0:
            continue
        # Jitter x positions
        rng_jitter = np.random.default_rng(42 + i)
        jitter = rng_jitter.uniform(-0.2, 0.2, n)
        x_positions.extend((np.ones(n) * (i + 1) + jitter).tolist())
        y_values.extend(by_cat[cat])
        colors_list.extend([CATEGORY_COLORS[cat]] * n)

    ax.scatter(x_positions, y_values, c=colors_list, alpha=0.7, s=40,
               edgecolors="white", linewidth=0.3)

    # Add horizontal lines for median
    for i, cat in enumerate(CATEGORY_ORDER):
        vals = by_cat[cat]
        if vals:
            median = np.median(vals)
            ax.hlines(y=median, xmin=i + 0.6, xmax=i + 1.4,
                      colors="black", linewidth=2, linestyles="-")

    ax.set_xticks(range(1, len(CATEGORY_ORDER) + 1))
    ax.set_xticklabels(CATEGORY_ORDER)
    ax.set_ylabel("BISER Overlap Fraction (%)")
    ax.set_title("Individual purge_dups Interval BISER Overlap Fractions")
    ax.set_ylim(-5, 105)

    # Add n annotations
    for i, cat in enumerate(CATEGORY_ORDER):
        n = len(by_cat[cat])
        ax.text(i + 1, 102, f"n={n}", ha="center", fontsize=8)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    png_path = FIGURES_DIR / "interval_overlap_fractions.png"
    pdf_path = FIGURES_DIR / "interval_overlap_fractions.pdf"
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"  Figure 4: {png_path}")
    return str(png_path)


# ── Figure 5: Chromosome-by-category heatmap ───────────────────────────────

def figure5_chromosome_heatmap(chrom_data: List[Dict]) -> str:
    """Heatmap/dot plot showing chromosome-specific overlap."""
    # Build matrix: rows=chromosomes, columns=categories
    # Only include chromosomes with at least one purge_dups interval
    chroms_with_data = set()
    for d in chrom_data:
        chroms_with_data.add(d["chromosome"])

    # Sort chromosomes: chr1-chr28, chrX, chrY, then seq* numerically
    def chrom_sort_key(c):
        if c.startswith("chr") and c[3:].isdigit():
            return (0, int(c[3:]))
        elif c == "chrX":
            return (0, 29)
        elif c == "chrY":
            return (0, 30)
        elif c.startswith("seq") and c[3:].isdigit():
            return (1, int(c[3:]))
        else:
            return (2, int(c[3:]) if c[3:].isdigit() else 999)

    sorted_chroms = sorted(chroms_with_data, key=chrom_sort_key)

    # Build data matrix
    matrix = np.zeros((len(sorted_chroms), len(CATEGORY_ORDER)))
    matrix.fill(np.nan)

    for d in chrom_data:
        cat = d["category"]
        chrom = d["chromosome"]
        if cat in CATEGORY_ORDER and chrom in sorted_chroms:
            row = sorted_chroms.index(chrom)
            col = CATEGORY_ORDER.index(cat)
            matrix[row, col] = float(d["pct_purge_bp_in_biser"])

    # Filter to chromosomes with at least 3 categories having data
    # Or show top 30 chromosomes by data presence
    n_data_per_chrom = np.sum(~np.isnan(matrix), axis=1)
    top_indices = np.argsort(n_data_per_chrom)[::-1][:35]
    matrix = matrix[top_indices]
    top_chroms = [sorted_chroms[i] for i in top_indices]

    fig, ax = plt.subplots(figsize=(8, max(6, len(top_chroms) * 0.3)))

    im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)

    ax.set_xticks(range(len(CATEGORY_ORDER)))
    ax.set_xticklabels(CATEGORY_ORDER, rotation=45, ha="right")
    ax.set_yticks(range(len(top_chroms)))
    ax.set_yticklabels(top_chroms, fontsize=7)

    # Add text annotations
    for i in range(len(top_chroms)):
        for j in range(len(CATEGORY_ORDER)):
            val = matrix[i, j]
            if not np.isnan(val):
                text_color = "white" if val > 60 else "black"
                ax.text(j, i, f"{val:.0f}%", ha="center", va="center",
                        fontsize=6, color=text_color)

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("BISER Overlap (%)")

    ax.set_title("Chromosome × Category BISER Overlap")
    fig.tight_layout()

    png_path = FIGURES_DIR / "chromosome_category_overlap.png"
    pdf_path = FIGURES_DIR / "chromosome_category_overlap.pdf"
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"  Figure 5: {png_path}")
    return str(png_path)


# ── Figure 6: Observed vs expected overlap ─────────────────────────────────

def figure6_observed_vs_expected(perm_results: List[Dict]) -> str:
    """Bar chart comparing observed vs expected (mean randomized) overlap."""
    fig, ax = plt.subplots(figsize=(7, 5))

    perm_by_cat = {r["category"]: r for r in perm_results}
    categories = [c for c in CATEGORY_ORDER if c in perm_by_cat]

    observed = [float(perm_by_cat[cat]["observed_overlap_mb"]) for cat in categories]
    expected = [float(perm_by_cat[cat]["null_mean_overlap_bp"]) / 1e6 for cat in categories]
    lower_ci = [float(perm_by_cat[cat]["null_lower_95_bp"]) / 1e6 for cat in categories]
    upper_ci = [float(perm_by_cat[cat]["null_upper_95_bp"]) / 1e6 for cat in categories]

    x = np.arange(len(categories))
    width = 0.35

    bars_obs = ax.bar(x - width/2, observed, width, color="darkred", alpha=0.8,
                       edgecolor="white", linewidth=0.5, label="Observed")
    bars_exp = ax.bar(x + width/2, expected, width, color="grey", alpha=0.6,
                       edgecolor="white", linewidth=0.5, label="Expected (null mean)")

    # Error bars on expected
    yerr_lower = [e - l for e, l in zip(expected, lower_ci)]
    yerr_upper = [u - e for e, u in zip(expected, upper_ci)]
    ax.errorbar(x + width/2, expected, yerr=[yerr_lower, yerr_upper],
                fmt="none", color="black", capsize=3, linewidth=1, label="95% CI")

    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_ylabel("BISER Overlap (Mb)")
    ax.set_title("Observed vs Expected BISER Overlap by Category")
    ax.legend(fontsize=8, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Significance annotations
    for i, cat in enumerate(categories):
        fdr = float(perm_by_cat[cat]["fdr_bh"])
        sig = ""
        if fdr < 0.05:
            sig = "*"
        if fdr < 0.01:
            sig = "**"
        if fdr < 0.001:
            sig = "***"
        if sig:
            max_h = max(observed[i], upper_ci[i])
            ax.text(i, max_h * 1.05, sig, ha="center", fontsize=14, fontweight="bold")

    fig.tight_layout()

    png_path = FIGURES_DIR / "observed_vs_expected_overlap.png"
    pdf_path = FIGURES_DIR / "observed_vs_expected_overlap.pdf"
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"  Figure 6: {png_path}")
    return str(png_path)


# ── Main ────────────────────────────────────────────────────────────────────

def main() -> None:
    """Generate all figures from cached results."""
    print("=" * 60)
    print("GENERATING FIGURES")
    print("=" * 60)

    # Load data
    obs = load_observed_table()
    perm_results = load_permutation_table()
    null_dists = load_null_distributions()
    chrom_data = load_chromosome_table()
    interval_details = load_interval_details()

    n_perm = perm_results[0]["n_permutations"] if perm_results else "?"
    run_mode = perm_results[0]["run_mode"] if perm_results else "?"
    print(f"  Run mode: {run_mode}, Permutations: {n_perm}")
    print(f"  Categories: {len(perm_results)}")

    # Generate figures
    figure1_observed_overlap(obs)
    figure2_fold_enrichment(perm_results)
    figure3_null_distributions(perm_results, null_dists)
    figure4_interval_overlap_fractions(interval_details)
    figure5_chromosome_heatmap(chrom_data)
    figure6_observed_vs_expected(perm_results)

    print("\n" + "=" * 60)
    print("ALL FIGURES GENERATED")
    print(f"Output directory: {FIGURES_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
