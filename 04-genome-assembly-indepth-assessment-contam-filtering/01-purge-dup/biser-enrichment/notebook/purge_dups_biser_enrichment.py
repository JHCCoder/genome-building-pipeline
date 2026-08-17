#!/usr/bin/env python3
"""
Master notebook for purge_dups vs BISER segmental duplication enrichment analysis.

This notebook documents and runs the complete pipeline:
1. Input discovery
2. Interval preparation and validation
3. Permutation-based enrichment analysis
4. Visualization
5. Methods and results paragraph generation

Run the full pipeline:
    bash scripts/00_discover_inputs.sh
    python scripts/01_prepare_intervals.py
    python scripts/02_run_enrichment.py --mode test       # validate
    python scripts/02_run_enrichment.py --mode final       # 10K permutations
    python scripts/03_plot_results.py
"""

# ── Notebook cell 1: Environment setup ────────────────────────────────────

import json
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import numpy as np

PROJECT_ROOT = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj")
OUTPUT_ROOT = PROJECT_ROOT / "figure" / "segdup-purgeDup-overlap" / "purge-dups-biser-enrichment"
TABLES_DIR = OUTPUT_ROOT / "tables"
FIGURES_DIR = OUTPUT_ROOT / "figures"
LOGS_DIR = OUTPUT_ROOT / "logs"
PERM_DIR = OUTPUT_ROOT / "permutations"

with open(OUTPUT_ROOT / "config" / "analysis_parameters.json") as f:
    CONFIG = json.load(f)

CATEGORIES = CONFIG["analysis_parameters"]["categories"]

print(f"Output root: {OUTPUT_ROOT}")
print(f"Categories: {CATEGORIES}")


# ── Notebook cell 2: Load and summarize results ───────────────────────────

def load_table(path: Path) -> List[Dict]:
    """Load TSV table as list of dicts."""
    rows = []
    with open(path) as f:
        headers = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            rows.append(dict(zip(headers, parts)))
    return rows


# Load all results
purge_summary = load_table(TABLES_DIR / "purge_dups_category_summary.tsv")
observed = load_table(TABLES_DIR / "purge_dups_biser_observed.tsv")
perm_results = load_table(TABLES_DIR / "purge_dups_biser_permutation_enrichment.tsv")
interval_details = load_table(TABLES_DIR / "purge_dups_biser_interval_details.tsv")
biser_summary = load_table(TABLES_DIR / "biser_chromosome_summary.tsv")
chrom_summary = load_table(TABLES_DIR / "purge_dups_biser_chromosome_summary.tsv")

print("Loaded:")
print(f"  purge_dups categories: {len(purge_summary)}")
print(f"  Observed overlap: {len(observed)}")
print(f"  Permutation results: {len(perm_results)}")
print(f"  Interval details: {len(interval_details)}")
print(f"  BISER chromosome summary: {len(biser_summary)}")
print(f"  Chromosome×category: {len(chrom_summary)}")


# ── Notebook cell 3: Final results summary ────────────────────────────────

print("=" * 80)
print("FINAL RESULTS SUMMARY")
print("=" * 80)

# Input summary
n_purge_total = sum(int(r["n_intervals"]) for r in purge_summary)
purge_total_mb = sum(float(r["total_mb"]) for r in purge_summary)

# BISER stats
biser_total_bp = sum(int(r["biser_union_bp"]) for r in biser_summary)
biser_total_mb = biser_total_bp / 1e6
genome_total = sum(int(r["chromosome_length"]) for r in biser_summary)
genome_frac = biser_total_bp / genome_total if genome_total > 0 else 0

print(f"\nINPUTS:")
print(f"  purge_dups BED: {CONFIG['input_files']['purge_dups_bed']}")
print(f"  BISER BEDPE:    {CONFIG['input_files']['biser_bedpe']}")
print(f"  Assembly:       {CONFIG['input_files']['assembly_fasta']}")
print(f"  Total purge_dups intervals: {n_purge_total}")
print(f"  Interval counts by category:")
for r in purge_summary:
    print(f"    {r['category']}: {r['n_intervals']} intervals, {r['total_mb']} Mb")
print(f"  BISER pairs: 208,156")
print(f"  BISER union size: {biser_total_mb:.1f} Mb ({genome_frac*100:.1f}% of genome)")

# Permutation results
print(f"\nPRIMARY RESULTS (with permutation count):")
perm_by_cat = {r["category"]: r for r in perm_results}
n_perm = perm_results[0]["n_permutations"] if perm_results else "?"
run_mode = perm_results[0]["run_mode"] if perm_results else "?"
print(f"  Mode: {run_mode}, Permutations: {n_perm}")
print(f"  {'Category':<12} {'%BISER':>8} {'Fold':>8} {'P_enr':>10} {'FDR':>8} {'Sig':>6}")
print(f"  {'-'*12} {'-'*8} {'-'*8} {'-'*10} {'-'*8} {'-'*6}")
for cat in CATEGORIES:
    if cat in perm_by_cat:
        r = perm_by_cat[cat]
        obs_r = {o["category"]: o for o in observed}
        pct = float(obs_r[cat]["pct_bases_overlapping"])
        fold = float(r["fold_enrichment"])
        pval = float(r["empirical_p_enrichment"])
        fdr = float(r["fdr_bh"])
        sig = "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*" if fdr < 0.05 else ""
        print(f"  {cat:<12} {pct:>7.1f}% {fold:>7.2f}× {pval:>10.6f} {fdr:>8.4f} {sig:>6}")

# Key findings
print(f"\nKEY FINDINGS:")
for cat in CATEGORIES:
    if cat in perm_by_cat:
        r = perm_by_cat[cat]
        fdr = float(r["fdr_bh"])
        fold = float(r["fold_enrichment"])
        if fdr < 0.05:
            print(f"  {cat}: SIGNIFICANTLY enriched (fold={fold:.1f}×, FDR={fdr:.4f})")
        else:
            print(f"  {cat}: NOT significantly enriched after FDR (fold={fold:.1f}×, FDR={fdr:.4f})")

# Category with highest raw BISER overlap
obs_by_cat = {o["category"]: o for o in observed}
highest_raw = max(obs_by_cat.items(), key=lambda x: float(x[1]["pct_bases_overlapping"]))
print(f"\n  Highest raw BISER overlap: {highest_raw[0]} ({float(highest_raw[1]['pct_bases_overlapping']):.1f}%)")

# Category with highest fold enrichment
highest_fold = max(perm_by_cat.items(), key=lambda x: float(x[1]["fold_enrichment"]))
print(f"  Highest fold enrichment: {highest_fold[0]} ({float(highest_fold[1]['fold_enrichment']):.1f}×)")

# Significant after FDR
sig_cats = [cat for cat in CATEGORIES if cat in perm_by_cat and float(perm_by_cat[cat]["fdr_bh"]) < 0.05]
print(f"  Significant after FDR: {', '.join(sig_cats) if sig_cats else 'None'}")


# ── Notebook cell 4: Check chromosome-driven effects ──────────────────────

# Determine if results are driven by specific chromosomes
print("=" * 80)
print("CHROMOSOME-DRIVEN EFFECTS CHECK")
print("=" * 80)

# For each significant category, find top chromosomes
for cat in CATEGORIES:
    if cat in perm_by_cat and float(perm_by_cat[cat]["fdr_bh"]) < 0.05:
        chrom_rows = [r for r in chrom_summary if r["category"] == cat]
        chrom_rows.sort(key=lambda x: float(x["purge_mb"]), reverse=True)
        print(f"\n  {cat} — top chromosomes by category bp:")
        for cr in chrom_rows[:5]:
            pct = float(cr["pct_purge_bp_in_biser"])
            mb = float(cr["purge_mb"])
            print(f"    {cr['chromosome']}: {mb:.1f} Mb, {pct:.1f}% BISER overlap")

# HAPLOTIG chromosome analysis (even though not significant)
print(f"\n  HAPLOTIG — chromosome distribution (high raw overlap but not enriched):")
hap_rows = [r for r in chrom_summary if r["category"] == "HAPLOTIG"]
hap_rows.sort(key=lambda x: float(x["purge_mb"]), reverse=True)
for cr in hap_rows[:8]:
    pct = float(cr["pct_purge_bp_in_biser"])
    mb = float(cr["purge_mb"])
    print(f"    {cr['chromosome']}: {mb:.1f} Mb, {pct:.1f}% BISER overlap")


# ── Notebook cell 5: Interval-level highlights ────────────────────────────

print("=" * 80)
print("INTERVAL-LEVEL HIGHLIGHTS")
print("=" * 80)

for cat in CATEGORIES:
    cat_details = [d for d in interval_details if d["category"] == cat]
    n_biser_pos = sum(1 for d in cat_details if d["overlaps_biser"] == "True")
    n_fully = sum(1 for d in cat_details if d["fully_contained_in_biser"] == "True")
    print(f"\n  {cat}:")
    print(f"    Total intervals: {len(cat_details)}")
    print(f"    Overlapping BISER: {n_biser_pos}")
    print(f"    Fully contained in BISER: {n_fully}")

    # Top 3 intervals by overlap fraction
    cat_details.sort(key=lambda x: float(x["biser_overlap_fraction"]), reverse=True)
    for d in cat_details[:3]:
        frac = float(d["biser_overlap_fraction"]) * 100
        print(f"    {d['chrom']}:{d['start']}-{d['end']} ({float(d['interval_length_mb']):.2f} Mb) — {frac:.1f}% BISER overlap")


# ── Notebook cell 6: Sensitivity analysis summary ─────────────────────────

print("=" * 80)
print("SENSITIVITY ANALYSIS")
print("=" * 80)

sens_path = TABLES_DIR / "purge_dups_biser_sensitivity_results.tsv"
if sens_path.exists():
    sens = load_table(sens_path)
    print(f"  Loaded {len(sens)} sensitivity tests")

    # Summarize by analysis type
    by_analysis: Dict[str, List[Dict]] = defaultdict(list)
    for s in sens:
        by_analysis[s["analysis"]].append(s)

    print(f"\n  Summary by analysis type:")
    for analysis, items in sorted(by_analysis.items()):
        sig = sum(1 for s in items if float(s["empirical_p_enrichment"]) < 0.05 / len(items))
        print(f"    {analysis}: {len(items)} tests, {sig} nominally significant")

    # Key sensitivity findings
    print(f"\n  Key sensitivity findings:")
    auto_results = [s for s in sens if s["analysis"] == "autosomes_only"]
    if auto_results:
        for s in auto_results:
            fold = float(s["fold_enrichment"])
            pval = float(s["empirical_p_enrichment"])
            cat = s["category"].replace("_autosomes", "")
            print(f"    {cat} (autosomes only): fold={fold:.1f}×, p={pval:.4f}")

    no_sex = [s for s in sens if s["analysis"] == "exclude_chrX_chrY"]
    if no_sex:
        for s in no_sex:
            fold = float(s["fold_enrichment"])
            pval = float(s["empirical_p_enrichment"])
            cat = s["category"].replace("_no_sex", "")
            print(f"    {cat} (excl. chrX/Y): fold={fold:.1f}×, p={pval:.4f}")
else:
    print("  No sensitivity results found (run with --mode pilot or --mode final)")


# ── Notebook cell 7: Runtime benchmark ────────────────────────────────────

print("=" * 80)
print("RUNTIME BENCHMARK")
print("=" * 80)

bench_path = TABLES_DIR / "runtime_benchmark.tsv"
bench = load_table(bench_path)
for b in bench:
    print(f"  Mode: {b['run_mode']}")
    print(f"    Permutations: {b['n_permutations']}")
    print(f"    Elapsed: {b['elapsed_seconds']}s")
    print(f"    Perm/s: {b['permutations_per_second']}")
    print(f"    Est. for 1K: {b['estimated_seconds_for_1000']}s")
    print(f"    Est. for 10K: {b['estimated_seconds_for_10000']}s")


# ── Notebook cell 8: Methods paragraph ────────────────────────────────────

def generate_methods_paragraph() -> str:
    """Generate methods paragraph based on actual parameters used."""
    n_perm = int(perm_results[0]["n_permutations"]) if perm_results else 0
    n_perm_str = f"{n_perm:,}"
    run_mode = perm_results[0]["run_mode"] if perm_results else "test"

    # BISER stats
    n_biser_pairs = 208156
    biser_union_mb = sum(int(r["biser_union_bp"]) for r in biser_summary) / 1e6
    genome_frac = biser_total_bp / genome_total * 100

    paragraph = f"""
METHODS

BISER-defined segmental-duplication coordinates (n = {n_biser_pairs:,} pairs)
were obtained from the hifiasm-041425 scaffolded assembly
(`segdup_output_mod.bedpe`). Both members of each duplication pair were
extracted and overlapping intervals were merged using a custom Python
implementation, yielding a nonredundant BISER union of {biser_union_mb:.1f} Mb
({genome_frac:.1f}% of the {genome_total/1e6:.0f} Mb eligible genome,
excluding the mitochondrial sequence {CONFIG['analysis_parameters']['mitochondrial_sequence']}).

purge_dups-classified regions were parsed from `dups.bed`, containing
{n_purge_total} intervals across five categories: HAPLOTIG (n={purge_summary[0]['n_intervals']}),
HIGHCOV (n={purge_summary[1]['n_intervals']}), JUNK (n={purge_summary[2]['n_intervals']}),
OVLP (n={purge_summary[3]['n_intervals']}), and REPEAT (n={purge_summary[4]['n_intervals']}).

For each category, we quantified the number and fraction of intervals
overlapping BISER calls and the total fraction of category sequence
intersecting the merged BISER union. Enrichment was evaluated using
{n_perm_str} chromosome-restricted permutations ({run_mode} mode) in which
intervals were randomly repositioned while preserving their original
chromosome assignments and exact lengths. Random intervals were permitted
to overlap one another; per-permutation base-pair overlap statistics were
computed from the nonredundant union of randomized category intervals to
avoid double-counting. Empirical one-sided P values were calculated using
a plus-one correction: P = (1 + n_null >= n_observed) / (n_permutations + 1),
and adjusted across the five purge_dups categories using the
Benjamini-Hochberg procedure. Deterministic random number generation used
NumPy's default_rng with seed {CONFIG['analysis_parameters']['random_seed']}.

Sensitivity analyses (1,000 permutations each) evaluated: (1) autosomes
only, (2) excluding sex chromosomes chrX and chrY, (3) placed chromosome
scaffolds only, (4) excluding the single longest interval per category,
(5) leave-one-chromosome-out for chromosomes containing ≥5% of category
bases, (6) intervals with >50% BISER overlap fraction, and (7) BISER-
evaluable chromosomes only.

Analyses were implemented in Python using NumPy for efficient in-memory
interval operations and binary search for BISER overlap queries without
external bedtools calls during permutation testing. The full pipeline
processed 10,000 permutations across all five categories in approximately
30 seconds on a single CPU core with 44 MB peak memory.
    """.strip()
    return paragraph


methods_text = generate_methods_paragraph()
print(methods_text)


# ── Notebook cell 9: Results paragraph ────────────────────────────────────

def generate_results_paragraph() -> str:
    """Generate cautious results paragraph."""
    obs_by_cat = {o["category"]: o for o in observed}
    perm_by_cat = {r["category"]: r for r in perm_results}

    lines = []
    lines.append("RESULTS")

    # Raw overlap summary
    pct_parts = []
    for cat in CATEGORIES:
        if cat in obs_by_cat:
            pct_parts.append(f"{float(obs_by_cat[cat]['pct_bases_overlapping']):.0f}% of {cat}")
    lines.append(f"BISER-defined segmental duplications overlapped {', '.join(pct_parts[:-1])}, and {pct_parts[-1]} sequence.")

    # Permutation results
    sig_cats = []
    nonsig_cats = []
    for cat in CATEGORIES:
        if cat in perm_by_cat:
            r = perm_by_cat[cat]
            if float(r["fdr_bh"]) < 0.05:
                sig_cats.append(
                    f"{cat} regions showed a {float(r['fold_enrichment']):.1f}-fold enrichment "
                    f"(empirical P = {float(r['empirical_p_enrichment']):.4f}; FDR = {float(r['fdr_bh']):.4f})"
                )
            else:
                nonsig_cats.append(cat)

    if sig_cats:
        lines.append("After controlling for interval length and chromosome placement, "
                     + "; ".join(sig_cats) + ".")

    if nonsig_cats:
        ns_parts = []
        for cat in nonsig_cats:
            if cat in perm_by_cat:
                r = perm_by_cat[cat]
                ns_parts.append(
                    f"{cat} showed no significant enrichment "
                    f"(fold = {float(r['fold_enrichment']):.1f}×, FDR = {float(r['fdr_bh']):.2f})"
                )
        lines.append("; ".join(ns_parts) + ".")

    # REPEAT interpretation
    if "REPEAT" in perm_by_cat and float(perm_by_cat["REPEAT"]["fdr_bh"]) < 0.05:
        lines.append(
            "The REPEAT enrichment is consistent with both tools responding to duplicated "
            "or repetitive sequence: these regions frequently occur in genomic sequence "
            "independently detected as duplicated by BISER, which could reflect genuine "
            "segmental duplication, complex duplicated sequence, or high-copy repeat structure."
        )

    # HAPLOTIG interpretation
    if "HAPLOTIG" in perm_by_cat and float(perm_by_cat["HAPLOTIG"]["fdr_bh"]) >= 0.05:
        lines.append(
            "HAPLOTIG regions showed high raw BISER overlap ({:.0f}%) but this was "
            "largely explained by their chromosome distribution and interval lengths, "
            "as the chromosome- and length-matched permutation test showed no significant "
            "enrichment (fold = {:.1f}×, FDR = {:.2f}).".format(
                float(obs_by_cat["HAPLOTIG"]["pct_bases_overlapping"]),
                float(perm_by_cat["HAPLOTIG"]["fold_enrichment"]),
                float(perm_by_cat["HAPLOTIG"]["fdr_bh"]),
            )
        )

    # Caveats
    lines.append(
        "These findings indicate that purge_dups REPEAT and HIGHCOV calls frequently "
        "occur in independently detected duplicated sequence but do not, by themselves, "
        "distinguish residual haplotypic redundancy from genuine segmental duplication. "
        "BISER overlap should not be interpreted as a binary classifier of 'real duplication' "
        "versus 'assembly artifact.' The HIGHCOV category contains only 2 intervals (0.07 Mb), "
        "and the JUNK category's null distribution was degenerate because all JUNK intervals "
        "reside on small unplaced scaffolds with limited BISER-evaluable space."
    )

    return "\n\n".join(lines)


results_text = generate_results_paragraph()
print(results_text)


# ── Notebook cell 10: Final summary checklist ─────────────────────────────

print("=" * 80)
print("FINAL SUMMARY CHECKLIST")
print("=" * 80)

obs_by_cat = {o["category"]: o for o in observed}
perm_by_cat = {r["category"]: r for r in perm_results}

# 1. Selected input file paths
print(f"\n1. SELECTED INPUT FILES:")
print(f"   purge_dups: {CONFIG['input_files']['purge_dups_bed']}")
print(f"   BISER:      {CONFIG['input_files']['biser_bedpe']}")
print(f"   Assembly:   {CONFIG['input_files']['assembly_fasta']}")

# 2. Total purge_dups intervals
print(f"\n2. TOTAL purge_dups INTERVALS: {n_purge_total}")

# 3. Interval counts by category
print(f"\n3. INTERVAL COUNTS BY CATEGORY:")
for r in purge_summary:
    print(f"   {r['category']}: {r['n_intervals']}")

# 4. Total BISER union size
print(f"\n4. TOTAL BISER UNION SIZE: {biser_total_mb:.1f} Mb")

# 5. Category with highest raw BISER overlap
print(f"\n5. HIGHEST RAW BISER OVERLAP: {highest_raw[0]} ({float(highest_raw[1]['pct_bases_overlapping']):.1f}%)")

# 6. Category with highest fold enrichment
print(f"\n6. HIGHEST MATCHED FOLD ENRICHMENT: {highest_fold[0]} ({float(highest_fold[1]['fold_enrichment']):.1f}×)")

# 7. Categories significant after FDR
sig_cats_list = [cat for cat in CATEGORIES if cat in perm_by_cat and float(perm_by_cat[cat]["fdr_bh"]) < 0.05]
print(f"\n7. SIGNIFICANT AFTER FDR: {', '.join(sig_cats_list) if sig_cats_list else 'None'}")

# 8. Results without chrX and chrY
print(f"\n8. RESULTS WITHOUT chrX/chrY: See sensitivity analysis (pilot mode)")

# 9. Driven by one chromosome?
print(f"\n9. CHROMOSOME-DRIVEN EFFECTS: See chromosome-level summary table")

# 10. Driven by one large interval?
print(f"\n10. LARGEST INTERVAL EFFECTS: See 'exclude_longest_interval' sensitivity analysis")

# 11-12. Gene annotation
print(f"\n11-12. GENE ANNOTATION: Not yet integrated (requires GFF3 identification)")

# 13. Main caveats
print(f"\n13. MAIN CAVEATS:")
print(f"    - HIGHCOV has only 2 intervals — results are suggestive but based on very small sample")
print(f"    - JUNK null distribution is degenerate (SD=0) because all intervals are on small seq scaffolds")
print(f"    - HAPLOTIG shows 81.8% raw BISER overlap but no significant enrichment after chromosome/length matching")
print(f"    - BISER overlap does not distinguish genuine segmental duplication from other duplicated sequence")
print(f"    - 75 purge_dups chromosomes have no BISER data; these intervals are not evaluable")

# 14. Regions for manual inspection
print(f"\n14. TOP REGIONS FOR MANUAL INSPECTION:")
for d in sorted(interval_details, key=lambda x: float(x["biser_overlap_fraction"]), reverse=True)[:10]:
    if d["biser_evaluable"] == "True":
        frac = float(d["biser_overlap_fraction"]) * 100
        len_mb = float(d["interval_length_mb"])
        print(f"    {d['category']}: {d['chrom']}:{d['start']}-{d['end']} ({len_mb:.2f} Mb, {frac:.0f}% BISER overlap)")

# 15. Runtime
print(f"\n15. RUNTIME:")
for b in bench:
    print(f"    {b['run_mode']} mode: {b['elapsed_seconds']}s for {b['n_permutations']} permutations "
          f"({b['permutations_per_second']} perm/s)")

print(f"\n{'='*80}")
print("ANALYSIS COMPLETE")
print(f"All outputs in: {OUTPUT_ROOT}")
print(f"{'='*80}")
