#!/usr/bin/env python3
"""
02_run_enrichment.py — Permutation-based enrichment analysis of purge_dups categories
for BISER-defined segmental duplications.

Implements efficient in-memory permutation testing with:
- Chromosome- and length-matched random interval placement
- Binary-search-based BISER overlap calculation
- Batching with progress reporting
- Deterministic random number generation
- Multiple run modes (test/pilot/final)
- Null distribution caching

Usage:
    python scripts/02_run_enrichment.py [--mode test|pilot|final]
"""

import argparse
import gzip
import json
import os
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj")
OUTPUT_ROOT = PROJECT_ROOT / "figure" / "segdup-purgeDup-overlap" / "purge-dups-biser-enrichment"
CONFIG_DIR = OUTPUT_ROOT / "config"
INTERMEDIATE_DIR = OUTPUT_ROOT / "intermediate"
TABLES_DIR = OUTPUT_ROOT / "tables"
PERM_DIR = OUTPUT_ROOT / "permutations"
LOGS_DIR = OUTPUT_ROOT / "logs"

with open(CONFIG_DIR / "analysis_parameters.json") as f:
    CONFIG = json.load(f)

CATEGORIES = CONFIG["analysis_parameters"]["categories"]
MITO_SEQ = CONFIG["analysis_parameters"]["mitochondrial_sequence"]
RANDOM_SEED = CONFIG["analysis_parameters"]["random_seed"]
BATCH_SIZE = CONFIG["analysis_parameters"]["batch_size"]
SENSITIVITY_N_PERM = CONFIG["analysis_parameters"]["sensitivity_n_permutations"]

PERMUTATIONS_BY_MODE = CONFIG["analysis_parameters"]["permutations_by_mode"]


# ── Data loading ───────────────────────────────────────────────────────────

def load_chromosome_sizes(fai_path: Path) -> Dict[str, int]:
    """Load chromosome sizes from FASTA index."""
    sizes: Dict[str, int] = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def load_bed(path: Path) -> List[Tuple[str, int, int]]:
    """Load BED3 file as list of (chrom, start, end)."""
    intervals: List[Tuple[str, int, int]] = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                intervals.append((parts[0], int(parts[1]), int(parts[2])))
    return intervals


def load_bed4(path: Path) -> List[Dict]:
    """Load BED4 file as list of dicts with chrom, start, end, name."""
    intervals: List[Dict] = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                intervals.append({
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "name": parts[3],
                })
    return intervals


def load_purge_dups_by_category() -> Dict[str, List[Dict]]:
    """Load purge_dups intervals grouped by category."""
    bed4_path = INTERMEDIATE_DIR / "purge_dups.categories.bed"
    intervals = load_bed4(bed4_path)
    by_cat: Dict[str, List[Dict]] = {c: [] for c in CATEGORIES}
    for iv in intervals:
        cat = iv["name"]
        if cat in CATEGORIES:
            by_cat[cat].append(iv)
    return by_cat


def load_biser_union_by_chrom() -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Load BISER union intervals as NumPy arrays per chromosome.

    Returns:
        Dict mapping chrom -> (starts_array, ends_array), both sorted.
    """
    union_path = INTERMEDIATE_DIR / "biser.segmental_duplication_union.bed"
    by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with open(union_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                by_chrom[parts[0]].append((int(parts[1]), int(parts[2])))

    result: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for chrom, ivs in by_chrom.items():
        ivs.sort(key=lambda x: x[0])
        starts = np.array([s for s, e in ivs], dtype=np.int64)
        ends = np.array([e for s, e in ivs], dtype=np.int64)
        result[chrom] = (starts, ends)
    return result


def load_eligible_genome() -> Dict[str, List[Tuple[int, int]]]:
    """Load eligible genome segments per chromosome.

    Returns:
        Dict mapping chrom -> list of (start, end) eligible segments.
    """
    eligible_path = INTERMEDIATE_DIR / "eligible_genome.bed"
    by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with open(eligible_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                by_chrom[parts[0]].append((int(parts[1]), int(parts[2])))
    return dict(by_chrom)


# ── Interval overlap calculation ───────────────────────────────────────────

def calculate_interval_overlap(
    chrom: str, start: int, end: int,
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
) -> int:
    """Calculate total bp overlap between one interval and BISER union.

    Uses binary search through sorted BISER intervals for efficiency.
    """
    if chrom not in biser_by_chrom:
        return 0
    b_starts, b_ends = biser_by_chrom[chrom]

    # Find index where b_starts < end (all BISER intervals starting before our end)
    right_idx = np.searchsorted(b_starts, end, side='left')
    # Find index where b_ends > start (all BISER intervals ending after our start)
    left_idx = np.searchsorted(b_ends, start, side='right')

    total = 0
    for i in range(left_idx, right_idx):
        overlap_start = max(start, int(b_starts[i]))
        overlap_end = min(end, int(b_ends[i]))
        if overlap_end > overlap_start:
            total += overlap_end - overlap_start
    return total


def calculate_union_overlap(
    intervals: List[Tuple[str, int, int]],
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
) -> int:
    """Calculate total nonredundant bp overlap between a set of intervals
    and the BISER union.

    First creates a nonredundant union of the input intervals per chromosome,
    then checks overlap against BISER.
    """
    # Build nonredundant union of input intervals per chromosome
    by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in intervals:
        by_chrom[chrom].append((start, end))

    total_overlap = 0
    for chrom, ivs in by_chrom.items():
        # Merge intervals
        ivs.sort(key=lambda x: x[0])
        merged: List[Tuple[int, int]] = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))

        # Calculate overlap with BISER
        for s, e in merged:
            total_overlap += calculate_interval_overlap(chrom, s, e, biser_by_chrom)
    return total_overlap


def count_overlapping_intervals(
    intervals: List[Tuple[str, int, int]],
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
) -> int:
    """Count how many intervals have at least 1 bp overlap with BISER."""
    count = 0
    for chrom, start, end in intervals:
        if calculate_interval_overlap(chrom, start, end, biser_by_chrom) > 0:
            count += 1
    return count


# ── Random interval placement ──────────────────────────────────────────────

def sample_valid_interval_starts(
    chrom_lengths: Dict[str, int],
    eligible_by_chrom: Dict[str, List[Tuple[int, int]]],
    intervals: List[Dict],
    rng: np.random.Generator,
) -> List[Tuple[str, int, int]]:
    """Randomly place intervals preserving chromosome and exact length.

    For each original interval, picks a uniform random valid start position
    on the same chromosome where the full interval fits within eligible space.

    Args:
        chrom_lengths: Dict mapping chromosome -> length in bp.
        eligible_by_chrom: Dict mapping chromosome -> list of (start, end) eligible segments.
        intervals: List of dicts with chrom, start, end, name keys.
        rng: NumPy random generator.

    Returns:
        List of (chrom, start, end) tuples for randomized intervals.
    """
    randomized: List[Tuple[str, int, int]] = []
    for iv in intervals:
        chrom = iv["chrom"]
        length = iv["end"] - iv["start"]

        if chrom not in eligible_by_chrom:
            # Chromosome not eligible — keep original position
            randomized.append((chrom, iv["start"], iv["end"]))
            continue

        eligible_segs = eligible_by_chrom[chrom]
        # Calculate valid start positions and their weights
        valid_ranges: List[Tuple[int, int, int]] = []  # (cumulative_weight, range_start, range_end)
        cumulative = 0
        for seg_start, seg_end in eligible_segs:
            max_start = seg_end - length
            if max_start > seg_start:
                n_positions = max_start - seg_start
                cumulative += n_positions
                valid_ranges.append((cumulative, seg_start, max_start))

        if not valid_ranges:
            # No valid placement — keep original
            randomized.append((chrom, iv["start"], iv["end"]))
            continue

        # Select a position uniformly from all valid positions
        rand_val = rng.integers(0, cumulative)
        for cum, range_start, range_end in valid_ranges:
            if rand_val < cum:
                # The selected position is within this range
                offset = rand_val - (cum - (range_end - range_start))
                new_start = range_start + offset
                randomized.append((chrom, new_start, new_start + length))
                break

    return randomized


# ── Observed overlap analysis ──────────────────────────────────────────────

def run_observed_analysis(
    category_intervals: Dict[str, List[Dict]],
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
    chrom_sizes: Dict[str, int],
    eligible_by_chrom: Dict[str, List[Tuple[int, int]]],
) -> Tuple[Dict[str, Dict], List[Dict]]:
    """Calculate observed overlap statistics for each category and interval."""
    print("\n" + "=" * 60)
    print("OBSERVED OVERLAP ANALYSIS")
    print("=" * 60)

    # Eligible genome BISER fraction
    biser_union_path = INTERMEDIATE_DIR / "biser.segmental_duplication_union.bed"
    total_eligible = sum(
        v for k, v in chrom_sizes.items()
        if k != MITO_SEQ and k in eligible_by_chrom
    )
    total_biser = 0
    with open(biser_union_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                chrom = parts[0]
                if chrom != MITO_SEQ:
                    total_biser += int(parts[2]) - int(parts[1])

    genome_biser_frac = total_biser / total_eligible if total_eligible > 0 else 0
    print(f"  Eligible genome BISER fraction: {genome_biser_frac:.4f} ({genome_biser_frac*100:.2f}%)")

    # Category-level observed analysis
    cat_results: Dict[str, Dict] = {}
    all_interval_details: List[Dict] = []

    for cat in CATEGORIES:
        ivs = category_intervals[cat]
        total_bp = sum(iv["end"] - iv["start"] for iv in ivs)

        # Nonredundant category union
        cat_union: List[Tuple[str, int, int]] = [
            (iv["chrom"], iv["start"], iv["end"]) for iv in ivs
        ]
        # Merge them
        by_chrom: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        for chrom, start, end in cat_union:
            by_chrom[chrom].append((start, end))
        cat_union_merged: List[Tuple[str, int, int]] = []
        for chrom, segs in by_chrom.items():
            segs.sort(key=lambda x: x[0])
            merged: List[Tuple[int, int]] = []
            for s, e in segs:
                if merged and s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            for s, e in merged:
                cat_union_merged.append((chrom, s, e))

        union_overlap = calculate_union_overlap(cat_union_merged, biser_by_chrom)
        n_overlapping = count_overlapping_intervals(cat_union, biser_by_chrom)
        n_fully_contained = 0
        n_partial = 0

        # Per-interval analysis
        for iv in ivs:
            chrom, start, end = iv["chrom"], iv["start"], iv["end"]
            length = end - start
            overlap = calculate_interval_overlap(chrom, start, end, biser_by_chrom)
            overlap_frac = overlap / length if length > 0 else 0

            # Overlap bin
            if overlap_frac == 0:
                overlap_bin = "0%"
            elif overlap_frac <= 0.10:
                overlap_bin = ">0-10%"
            elif overlap_frac <= 0.25:
                overlap_bin = ">10-25%"
            elif overlap_frac <= 0.50:
                overlap_bin = ">25-50%"
            elif overlap_frac <= 0.75:
                overlap_bin = ">50-75%"
            elif overlap_frac < 1.0:
                overlap_bin = ">75-99%"
            else:
                overlap_bin = "100%"

            biser_evaluable = chrom in biser_by_chrom

            if overlap == length and length > 0:
                n_fully_contained += 1
            elif overlap > 0:
                n_partial += 1

            detail = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "category": cat,
                "interval_length_bp": length,
                "interval_length_mb": length / 1e6,
                "biser_overlap_bp": overlap,
                "biser_overlap_fraction": overlap_frac,
                "overlaps_biser": overlap > 0,
                "fully_contained_in_biser": overlap == length and length > 0,
                "overlap_bin": overlap_bin,
                "biser_evaluable": biser_evaluable,
            }
            all_interval_details.append(detail)

        cat_overlap_frac = union_overlap / total_bp if total_bp > 0 else 0
        naive_fe = cat_overlap_frac / genome_biser_frac if genome_biser_frac > 0 else float('inf')

        cat_results[cat] = {
            "n_intervals": len(ivs),
            "total_bp": total_bp,
            "total_mb": total_bp / 1e6,
            "n_overlapping": n_overlapping,
            "pct_overlapping": n_overlapping / len(ivs) * 100 if ivs else 0,
            "n_fully_contained": n_fully_contained,
            "n_partial": n_partial,
            "biser_overlap_bp": union_overlap,
            "biser_overlap_mb": union_overlap / 1e6,
            "pct_bases_overlapping": cat_overlap_frac * 100,
            "non_biser_bp": total_bp - union_overlap,
            "naive_fold_enrichment": naive_fe,
            "eligible_genome_biser_frac": genome_biser_frac,
        }

        print(f"\n  {cat}:")
        print(f"    Intervals: {len(ivs)}")
        print(f"    Total bases: {total_bp:,} bp ({total_bp/1e6:.2f} Mb)")
        print(f"    Overlapping BISER: {n_overlapping}/{len(ivs)} ({n_overlapping/len(ivs)*100:.1f}%)")
        print(f"    BISER overlap bp: {union_overlap:,} ({cat_overlap_frac*100:.1f}%)")
        print(f"    Naive fold enrichment: {naive_fe:.2f}")

    # Write observed table
    obs_path = TABLES_DIR / "purge_dups_biser_observed.tsv"
    with open(obs_path, "w") as f:
        headers = [
            "category", "n_intervals", "total_bp", "total_mb",
            "n_overlapping", "pct_overlapping", "n_fully_contained",
            "n_partial", "biser_overlap_bp", "biser_overlap_mb",
            "pct_bases_overlapping", "non_biser_bp",
            "naive_fold_enrichment", "eligible_genome_biser_frac",
        ]
        f.write("\t".join(headers) + "\n")
        for cat in CATEGORIES:
            r = cat_results[cat]
            f.write("\t".join(str(r[h]) for h in headers[1:]) + "\n")
            # Prepend category
    # Rewrite with proper category column
    lines = []
    with open(obs_path) as f:
        lines = f.readlines()
    with open(obs_path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for cat in CATEGORIES:
            r = cat_results[cat]
            vals = [cat] + [str(r[h]) for h in headers[1:]]
            f.write("\t".join(vals) + "\n")
    print(f"\n  Wrote: {obs_path}")

    # Write interval details
    detail_path = TABLES_DIR / "purge_dups_biser_interval_details.tsv"
    detail_headers = [
        "chrom", "start", "end", "category", "interval_length_bp",
        "interval_length_mb", "biser_overlap_bp", "biser_overlap_fraction",
        "overlaps_biser", "fully_contained_in_biser", "overlap_bin",
        "biser_evaluable",
    ]
    with open(detail_path, "w") as f:
        f.write("\t".join(detail_headers) + "\n")
        for d in sorted(all_interval_details, key=lambda x: (CATEGORIES.index(x["category"]), x["chrom"], x["start"])):
            f.write("\t".join(str(d[h]) for h in detail_headers) + "\n")
    print(f"  Wrote: {detail_path}")

    return cat_results, all_interval_details


# ── Permutation test ───────────────────────────────────────────────────────

def run_permutation_batch(
    category_intervals: List[Dict],
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
    chrom_sizes: Dict[str, int],
    eligible_by_chrom: Dict[str, List[Tuple[int, int]]],
    n_perm: int,
    rng: np.random.Generator,
    batch_id: int = 0,
) -> Tuple[List[int], List[int]]:
    """Run a batch of permutations for one category.

    Args:
        category_intervals: List of interval dicts for this category.
        biser_by_chrom: BISER union intervals per chromosome.
        chrom_sizes: Chromosome sizes.
        eligible_by_chrom: Eligible segments per chromosome.
        n_perm: Number of permutations in this batch.
        rng: NumPy random generator.
        batch_id: Batch identifier for progress reporting.

    Returns:
        (overlap_bp_list, n_overlapping_list) — one value per permutation.
    """
    overlap_results: List[int] = []
    n_overlap_results: List[int] = []

    for i in range(n_perm):
        randomized = sample_valid_interval_starts(
            chrom_sizes, eligible_by_chrom, category_intervals, rng,
        )
        overlap = calculate_union_overlap(randomized, biser_by_chrom)
        n_ov = count_overlapping_intervals(randomized, biser_by_chrom)
        overlap_results.append(overlap)
        n_overlap_results.append(n_ov)

    return overlap_results, n_overlap_results


def run_permutation_test(
    category_name: str,
    category_intervals: List[Dict],
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
    chrom_sizes: Dict[str, int],
    eligible_by_chrom: Dict[str, List[Tuple[int, int]]],
    n_permutations: int,
    batch_size: int,
    rng: np.random.Generator,
    observed_overlap_bp: int,
    observed_n_overlapping: int,
) -> Dict:
    """Run full permutation test for one category.

    Returns dict with null distribution and statistics.
    """
    n_batches = (n_permutations + batch_size - 1) // batch_size

    all_overlaps: List[int] = []
    all_n_overlapping: List[int] = []

    t0 = time.time()
    for batch_i in range(n_batches):
        this_batch_size = min(batch_size, n_permutations - len(all_overlaps))
        overlaps, n_ovs = run_permutation_batch(
            category_intervals, biser_by_chrom, chrom_sizes,
            eligible_by_chrom, this_batch_size, rng, batch_i,
        )
        all_overlaps.extend(overlaps)
        all_n_overlapping.extend(n_ovs)

        elapsed = time.time() - t0
        done = len(all_overlaps)
        rate = done / elapsed if elapsed > 0 else 0
        eta = (n_permutations - done) / rate if rate > 0 else 0
        print(f"  {category_name} batch {batch_i+1}/{n_batches}: "
              f"{done}/{n_permutations} done, {rate:.0f} perm/s, ETA {eta:.0f}s")

    null_overlaps = np.array(all_overlaps, dtype=np.int64)
    null_n_ov = np.array(all_n_overlapping, dtype=np.int64)

    # Statistics for overlap bp
    null_mean = float(np.mean(null_overlaps))
    null_median = float(np.median(null_overlaps))
    null_sd = float(np.std(null_overlaps))
    null_lower = float(np.percentile(null_overlaps, 2.5))
    null_upper = float(np.percentile(null_overlaps, 97.5))

    # Empirical p-values (plus-one correction)
    n_gte = int(np.sum(null_overlaps >= observed_overlap_bp))
    n_lte = int(np.sum(null_overlaps <= observed_overlap_bp))
    p_enrichment = (1 + n_gte) / (n_permutations + 1)
    p_depletion = (1 + n_lte) / (n_permutations + 1)

    # Fold enrichment
    fold_enrichment = observed_overlap_bp / null_mean if null_mean > 0 else float('inf')

    # z-score
    z_score = (observed_overlap_bp - null_mean) / null_sd if null_sd > 0 else float('nan')

    # Statistics for n_overlapping
    null_n_ov_mean = float(np.mean(null_n_ov))
    null_n_ov_sd = float(np.std(null_n_ov))

    return {
        "category": category_name,
        "observed_overlap_bp": observed_overlap_bp,
        "observed_overlap_mb": observed_overlap_bp / 1e6,
        "observed_n_overlapping": observed_n_overlapping,
        "null_mean_overlap_bp": null_mean,
        "null_median_overlap_bp": null_median,
        "null_sd_overlap_bp": null_sd,
        "null_lower_95_bp": null_lower,
        "null_upper_95_bp": null_upper,
        "fold_enrichment": fold_enrichment,
        "z_score": z_score,
        "empirical_p_enrichment": p_enrichment,
        "empirical_p_depletion": p_depletion,
        "null_n_ov_mean": null_n_ov_mean,
        "null_n_ov_sd": null_n_ov_sd,
        "null_overlaps": null_overlaps,
        "null_n_overlapping": null_n_ov,
        "n_permutations": n_permutations,
    }


# ── FDR correction ─────────────────────────────────────────────────────────

def benjamini_hochberg(pvalues: List[Tuple[str, float]]) -> List[Tuple[str, float, float]]:
    """Apply Benjamini-Hochberg FDR correction.

    Args:
        pvalues: List of (name, p_value) tuples.

    Returns:
        List of (name, p_value, fdr) tuples sorted by p-value.
    """
    sorted_pv = sorted(pvalues, key=lambda x: x[1])
    n = len(sorted_pv)
    fdr_values: List[float] = []
    for i, (name, pv) in enumerate(sorted_pv):
        # BH: p_i * n / rank_i, with cumulative min from the end
        fdr_values.append(min(pv * n / (i + 1), 1.0))

    # Ensure monotonicity (cumulative min from end)
    for i in range(n - 2, -1, -1):
        fdr_values[i] = min(fdr_values[i], fdr_values[i + 1])

    return [(name, pv, fdr) for (name, pv), fdr in zip(sorted_pv, fdr_values)]


# ── Sensitivity analyses ───────────────────────────────────────────────────

def run_sensitivity_analyses(
    category_intervals: Dict[str, List[Dict]],
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
    chrom_sizes: Dict[str, int],
    eligible_by_chrom: Dict[str, List[Tuple[int, int]]],
    n_permutations: int,
    batch_size: int,
    base_rng: np.random.Generator,
) -> List[Dict]:
    """Run runtime-conscious sensitivity analyses."""
    print("\n" + "=" * 60)
    print("SENSITIVITY ANALYSES")
    print(f"  Using {n_permutations} permutations per test")
    print("=" * 60)

    sensitivity_results: List[Dict] = []

    # 1. Autosomes only (exclude chrX, chrY, seq*)
    print("\n  [1/7] Autosomes only (exclude chrX, chrY, seq*) ...")
    auto_intervals = {}
    for cat in CATEGORIES:
        auto_intervals[cat] = [
            iv for iv in category_intervals[cat]
            if iv["chrom"].startswith("chr") and iv["chrom"] not in ("chrX", "chrY")
            and iv["chrom"] in eligible_by_chrom
        ]
    auto_chroms = {
        k: v for k, v in chrom_sizes.items()
        if k.startswith("chr") and k not in ("chrX", "chrY")
    }
    auto_eligible = {
        k: v for k, v in eligible_by_chrom.items()
        if k.startswith("chr") and k not in ("chrX", "chrY")
    }
    for cat in CATEGORIES:
        if len(auto_intervals[cat]) == 0:
            continue
        ivs = auto_intervals[cat]
        cat_union = [(iv["chrom"], iv["start"], iv["end"]) for iv in ivs]
        obs_ov = calculate_union_overlap(cat_union, biser_by_chrom)
        obs_n = count_overlapping_intervals(cat_union, biser_by_chrom)
        seed = base_rng.integers(0, 2**31 - 1)
        rng = np.random.default_rng(seed)
        result = run_permutation_test(
            f"{cat}_autosomes", ivs, biser_by_chrom, auto_chroms,
            auto_eligible, n_permutations, batch_size, rng, obs_ov, obs_n,
        )
        result["analysis"] = "autosomes_only"
        sensitivity_results.append(result)

    # 2. Excluding chrX and chrY
    print("\n  [2/7] Excluding chrX and chrY ...")
    no_sex_intervals = {}
    for cat in CATEGORIES:
        no_sex_intervals[cat] = [
            iv for iv in category_intervals[cat]
            if iv["chrom"] not in ("chrX", "chrY") and iv["chrom"] in eligible_by_chrom
        ]
    no_sex_chroms = {k: v for k, v in chrom_sizes.items() if k not in ("chrX", "chrY")}
    no_sex_eligible = {k: v for k, v in eligible_by_chrom.items() if k not in ("chrX", "chrY")}
    for cat in CATEGORIES:
        if len(no_sex_intervals[cat]) == 0:
            continue
        ivs = no_sex_intervals[cat]
        cat_union = [(iv["chrom"], iv["start"], iv["end"]) for iv in ivs]
        obs_ov = calculate_union_overlap(cat_union, biser_by_chrom)
        obs_n = count_overlapping_intervals(cat_union, biser_by_chrom)
        seed = base_rng.integers(0, 2**31 - 1)
        rng = np.random.default_rng(seed)
        result = run_permutation_test(
            f"{cat}_no_sex", ivs, biser_by_chrom, no_sex_chroms,
            no_sex_eligible, n_permutations, batch_size, rng, obs_ov, obs_n,
        )
        result["analysis"] = "exclude_chrX_chrY"
        sensitivity_results.append(result)

    # 3. Excluding unplaced scaffolds
    print("\n  [3/7] Excluding unplaced scaffolds ...")
    placed_intervals = {}
    for cat in CATEGORIES:
        placed_intervals[cat] = [
            iv for iv in category_intervals[cat]
            if iv["chrom"].startswith("chr") and iv["chrom"] in eligible_by_chrom
        ]
    placed_chroms = {k: v for k, v in chrom_sizes.items() if k.startswith("chr")}
    placed_eligible = {k: v for k, v in eligible_by_chrom.items() if k.startswith("chr")}
    for cat in CATEGORIES:
        if len(placed_intervals[cat]) == 0:
            continue
        ivs = placed_intervals[cat]
        cat_union = [(iv["chrom"], iv["start"], iv["end"]) for iv in ivs]
        obs_ov = calculate_union_overlap(cat_union, biser_by_chrom)
        obs_n = count_overlapping_intervals(cat_union, biser_by_chrom)
        seed = base_rng.integers(0, 2**31 - 1)
        rng = np.random.default_rng(seed)
        result = run_permutation_test(
            f"{cat}_placed_only", ivs, biser_by_chrom, placed_chroms,
            placed_eligible, n_permutations, batch_size, rng, obs_ov, obs_n,
        )
        result["analysis"] = "placed_scaffolds_only"
        sensitivity_results.append(result)

    # 4. Excluding the single longest interval per category
    print("\n  [4/7] Excluding longest interval per category ...")
    for cat in CATEGORIES:
        ivs = sorted(category_intervals[cat], key=lambda x: x["end"] - x["start"], reverse=True)
        if len(ivs) <= 1:
            continue
        trimmed = ivs[1:]  # Remove longest
        trimmed = [iv for iv in trimmed if iv["chrom"] in eligible_by_chrom]
        if len(trimmed) == 0:
            continue
        cat_union = [(iv["chrom"], iv["start"], iv["end"]) for iv in trimmed]
        obs_ov = calculate_union_overlap(cat_union, biser_by_chrom)
        obs_n = count_overlapping_intervals(cat_union, biser_by_chrom)
        seed = base_rng.integers(0, 2**31 - 1)
        rng = np.random.default_rng(seed)
        result = run_permutation_test(
            f"{cat}_no_longest", trimmed, biser_by_chrom, chrom_sizes,
            eligible_by_chrom, n_permutations, batch_size, rng, obs_ov, obs_n,
        )
        result["analysis"] = "exclude_longest_interval"
        sensitivity_results.append(result)

    # 5. Leave-one-chromosome-out for top categories
    # Only test chromosomes with >= 5% of the category's total bases
    print("\n  [5/7] Leave-one-chromosome-out (chromosomes with >=5% of category bp) ...")
    for cat in CATEGORIES:
        ivs = category_intervals[cat]
        total_cat_bp = sum(iv["end"] - iv["start"] for iv in ivs)
        chrom_bp: Dict[str, int] = defaultdict(int)
        for iv in ivs:
            chrom_bp[iv["chrom"]] += iv["end"] - iv["start"]
        major_chroms = [c for c, bp in chrom_bp.items() if bp >= 0.05 * total_cat_bp]
        print(f"    {cat}: testing {len(major_chroms)}/{len(chrom_bp)} chromosomes")
        for chrom in sorted(major_chroms):
            loo_ivs = [iv for iv in ivs if iv["chrom"] != chrom and iv["chrom"] in eligible_by_chrom]
            if len(loo_ivs) == 0:
                continue
            cat_union = [(iv["chrom"], iv["start"], iv["end"]) for iv in loo_ivs]
            obs_ov = calculate_union_overlap(cat_union, biser_by_chrom)
            obs_n = count_overlapping_intervals(cat_union, biser_by_chrom)
            seed = base_rng.integers(0, 2**31 - 1)
            rng = np.random.default_rng(seed)
            result = run_permutation_test(
                f"{cat}_without_{chrom}", loo_ivs, biser_by_chrom, chrom_sizes,
                eligible_by_chrom, n_permutations, batch_size, rng, obs_ov, obs_n,
            )
            result["analysis"] = f"leave_one_out_{chrom}"
            sensitivity_results.append(result)

    # 6. Reciprocal threshold analysis
    print("\n  [6/7] High-BISER-overlap intervals ...")
    for cat in CATEGORIES:
        high_ov_ivs = [
            iv for iv in category_intervals[cat]
            if iv["chrom"] in eligible_by_chrom and
            calculate_interval_overlap(iv["chrom"], iv["start"], iv["end"], biser_by_chrom) /
            (iv["end"] - iv["start"]) > 0.5
        ]
        if len(high_ov_ivs) < 2:
            continue
        cat_union = [(iv["chrom"], iv["start"], iv["end"]) for iv in high_ov_ivs]
        obs_ov = calculate_union_overlap(cat_union, biser_by_chrom)
        obs_n = count_overlapping_intervals(cat_union, biser_by_chrom)
        seed = base_rng.integers(0, 2**31 - 1)
        rng = np.random.default_rng(seed)
        result = run_permutation_test(
            f"{cat}_high_biser_overlap", high_ov_ivs, biser_by_chrom, chrom_sizes,
            eligible_by_chrom, n_permutations, batch_size, rng, obs_ov, obs_n,
        )
        result["analysis"] = "high_biser_overlap_intervals"
        sensitivity_results.append(result)

    # 7. BISER-evaluable chromosomes only
    print("\n  [7/7] BISER-evaluable chromosomes only ...")
    biser_chroms = set(biser_by_chrom.keys())
    for cat in CATEGORIES:
        eval_ivs = [
            iv for iv in category_intervals[cat]
            if iv["chrom"] in biser_chroms and iv["chrom"] in eligible_by_chrom
        ]
        if len(eval_ivs) == 0:
            continue
        cat_union = [(iv["chrom"], iv["start"], iv["end"]) for iv in eval_ivs]
        obs_ov = calculate_union_overlap(cat_union, biser_by_chrom)
        obs_n = count_overlapping_intervals(cat_union, biser_by_chrom)
        eval_chroms = {k: v for k, v in chrom_sizes.items() if k in biser_chroms}
        eval_eligible = {k: v for k, v in eligible_by_chrom.items() if k in biser_chroms}
        seed = base_rng.integers(0, 2**31 - 1)
        rng = np.random.default_rng(seed)
        result = run_permutation_test(
            f"{cat}_biser_evaluable", eval_ivs, biser_by_chrom, eval_chroms,
            eval_eligible, n_permutations, batch_size, rng, obs_ov, obs_n,
        )
        result["analysis"] = "biser_evaluable_only"
        sensitivity_results.append(result)

    return sensitivity_results


# ── Chromosome-level analysis ──────────────────────────────────────────────

def run_chromosome_analysis(
    category_intervals: Dict[str, List[Dict]],
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
    chrom_sizes: Dict[str, int],
) -> None:
    """Generate chromosome-by-category summary."""
    print("\n" + "=" * 60)
    print("CHROMOSOME-LEVEL ANALYSIS")
    print("=" * 60)

    chr_summary_path = TABLES_DIR / "purge_dups_biser_chromosome_summary.tsv"
    with open(chr_summary_path, "w") as f:
        f.write("chromosome\tcategory\tn_purge_regions\tpurge_bp\tpurge_mb\t"
                "biser_overlap_bp\tbiser_overlap_mb\tpct_purge_bp_in_biser\t"
                "chromosome_biser_fraction\tchromosome_length\n")

        for chrom in sorted(chrom_sizes.keys(), key=lambda x: (
            0 if x.startswith("chr") and x[3:].isdigit() else
            1 if x in ("chrX", "chrY") else 2,
            int(x[3:]) if x.startswith("chr") and x[3:].isdigit() else
            (29 if x == "chrX" else 30 if x == "chrY" else int(x[3:]) if x.startswith("seq") and x[3:].isdigit() else 999)
        )):
            chrom_len = chrom_sizes[chrom]
            # Chromosome BISER fraction
            chrom_biser_bp = 0
            if chrom in biser_by_chrom:
                starts, ends = biser_by_chrom[chrom]
                chrom_biser_bp = int(np.sum(ends - starts))
            chrom_biser_frac = chrom_biser_bp / chrom_len if chrom_len > 0 else 0

            for cat in CATEGORIES:
                ivs = [iv for iv in category_intervals[cat] if iv["chrom"] == chrom]
                if not ivs:
                    continue
                total_bp = sum(iv["end"] - iv["start"] for iv in ivs)
                # Union overlap
                cat_union = [(iv["chrom"], iv["start"], iv["end"]) for iv in ivs]
                overlap = calculate_union_overlap(cat_union, biser_by_chrom)
                pct = overlap / total_bp * 100 if total_bp > 0 else 0

                f.write(
                    f"{chrom}\t{cat}\t{len(ivs)}\t{total_bp}\t{total_bp/1e6:.2f}\t"
                    f"{overlap}\t{overlap/1e6:.2f}\t{pct:.2f}\t"
                    f"{chrom_biser_frac:.4f}\t{chrom_len}\n"
                )

    print(f"  Wrote: {chr_summary_path}")


# ── Runtime benchmarking ───────────────────────────────────────────────────

def benchmark_permutations(
    category_intervals: Dict[str, List[Dict]],
    biser_by_chrom: Dict[str, Tuple[np.ndarray, np.ndarray]],
    chrom_sizes: Dict[str, int],
    eligible_by_chrom: Dict[str, List[Tuple[int, int]]],
    rng: np.random.Generator,
) -> Dict:
    """Benchmark permutation speed and estimate runtime for larger runs."""
    print("\n" + "=" * 60)
    print("RUNTIME BENCHMARK")
    print("=" * 60)

    # Pick the category with the most intervals for worst-case benchmark
    largest_cat = max(CATEGORIES, key=lambda c: len(category_intervals[c]))
    ivs = category_intervals[largest_cat]
    n_bench = 50

    t0 = time.time()
    overlaps, n_ovs = run_permutation_batch(
        ivs, biser_by_chrom, chrom_sizes, eligible_by_chrom, n_bench, rng, 0,
    )
    elapsed = time.time() - t0
    rate = n_bench / elapsed

    estimates = {}
    for mode, n in PERMUTATIONS_BY_MODE.items():
        est = n / rate
        estimates[mode] = est
        print(f"  {largest_cat} ({len(ivs)} intervals): {rate:.0f} perm/s")
        print(f"    Estimated for {n:,} perms: {est:.1f}s ({est/60:.1f} min)")
        print(f"    Estimated for all 5 categories: {est*5:.1f}s ({est*5/60:.1f} min)")

    return {
        "benchmark_category": largest_cat,
        "n_intervals": len(ivs),
        "n_benchmark_perms": n_bench,
        "elapsed_seconds": elapsed,
        "permutations_per_second": rate,
        "estimates": estimates,
    }


# ── Main ────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Run purge_dups vs BISER enrichment analysis")
    parser.add_argument("--mode", choices=["test", "pilot", "final"], default="test",
                        help="Run mode: test (100), pilot (1000), final (10000)")
    parser.add_argument("--no-sensitivity", action="store_true",
                        help="Skip sensitivity analyses")
    args = parser.parse_args()

    mode = args.mode
    n_permutations = PERMUTATIONS_BY_MODE[mode]
    batch_size = BATCH_SIZE

    print("=" * 60)
    print("PURGE_DUPS vs BISER: Permutation Enrichment Analysis")
    print("=" * 60)
    print(f"  Mode: {mode}")
    print(f"  Permutations: {n_permutations:,}")
    print(f"  Batch size: {batch_size}")
    print(f"  Random seed: {RANDOM_SEED}")

    # Load data
    chrom_sizes = load_chromosome_sizes(CONFIG["input_files"]["assembly_fai"])
    category_intervals = load_purge_dups_by_category()
    biser_by_chrom = load_biser_union_by_chrom()
    eligible_by_chrom = load_eligible_genome()

    # Print summary
    total_intervals = sum(len(v) for v in category_intervals.values())
    print(f"\n  Loaded {total_intervals} purge_dups intervals across {len(CATEGORIES)} categories")
    for cat in CATEGORIES:
        print(f"    {cat}: {len(category_intervals[cat])} intervals")
    print(f"  BISER union: {sum(len(v[0]) for v in biser_by_chrom.values())} intervals "
          f"across {len(biser_by_chrom)} chromosomes")

    # Create RNG
    rng = np.random.default_rng(RANDOM_SEED)

    # Step 1: Observed overlap analysis
    cat_results, all_interval_details = run_observed_analysis(
        category_intervals, biser_by_chrom, chrom_sizes, eligible_by_chrom,
    )

    # Step 2: Benchmark
    bench = benchmark_permutations(
        category_intervals, biser_by_chrom, chrom_sizes, eligible_by_chrom, rng,
    )

    # Step 3: Primary permutation test
    print("\n" + "=" * 60)
    print(f"PRIMARY PERMUTATION TEST ({mode}: {n_permutations:,} permutations)")
    print("=" * 60)

    # Reset RNG for reproducibility
    rng = np.random.default_rng(RANDOM_SEED)

    all_results: List[Dict] = []
    for cat in CATEGORIES:
        ivs = category_intervals[cat]
        if not ivs:
            print(f"\n  {cat}: No intervals, skipping")
            continue
        print(f"\n  Running {cat} ({len(ivs)} intervals)...")
        obs = cat_results[cat]
        result = run_permutation_test(
            cat, ivs, biser_by_chrom, chrom_sizes, eligible_by_chrom,
            n_permutations, batch_size, rng,
            obs["biser_overlap_bp"], obs["n_overlapping"],
        )
        all_results.append(result)

    # Apply FDR correction
    pvalues = [(r["category"], r["empirical_p_enrichment"]) for r in all_results]
    fdr_corrected = benjamini_hochberg(pvalues)
    fdr_map = {name: fdr for name, pv, fdr in fdr_corrected}

    for r in all_results:
        r["fdr_bh"] = fdr_map.get(r["category"], 1.0)
        r["n_permutations"] = n_permutations
        r["run_mode"] = mode

    # Write permutation enrichment table
    enrich_path = TABLES_DIR / "purge_dups_biser_permutation_enrichment.tsv"
    enrich_headers = [
        "category", "observed_overlap_bp", "observed_overlap_mb",
        "observed_n_overlapping", "null_mean_overlap_bp", "null_median_overlap_bp",
        "null_sd_overlap_bp", "null_lower_95_bp", "null_upper_95_bp",
        "fold_enrichment", "z_score", "empirical_p_enrichment",
        "empirical_p_depletion", "fdr_bh", "n_permutations", "run_mode",
    ]
    with open(enrich_path, "w") as f:
        f.write("\t".join(enrich_headers) + "\n")
        for r in all_results:
            f.write("\t".join(str(r.get(h, "")) for h in enrich_headers) + "\n")
    print(f"\n  Wrote: {enrich_path}")

    # Print summary
    print(f"\n{'Category':<12} {'Obs(Mb)':>8} {'NullMn':>8} {'Null95Lo':>8} "
          f"{'Null95Hi':>8} {'Fold':>8} {'P(>)=)':>10} {'FDR':>8}")
    print("-" * 80)
    for r in all_results:
        print(f"{r['category']:<12} {r['observed_overlap_mb']:>7.2f} "
              f"{r['null_mean_overlap_bp']/1e6:>7.2f} "
              f"{r['null_lower_95_bp']/1e6:>7.2f} "
              f"{r['null_upper_95_bp']/1e6:>7.2f} "
              f"{r['fold_enrichment']:>7.2f} "
              f"{r['empirical_p_enrichment']:>10.6f} "
              f"{r['fdr_bh']:>8.4f}")

    # Save null distributions
    print(f"\n  Saving null distributions...")
    npz_path = PERM_DIR / "primary_null_distributions.npz"
    save_dict = {}
    for r in all_results:
        save_dict[f"{r['category']}_overlap"] = r["null_overlaps"]
        save_dict[f"{r['category']}_n_overlapping"] = r["null_n_overlapping"]
    np.savez_compressed(npz_path, **save_dict)
    print(f"  Wrote: {npz_path}")

    # Step 4: Run sensitivity analyses (only for pilot and final modes)
    if mode in ("pilot", "final") and not args.no_sensitivity:
        sens_n = SENSITIVITY_N_PERM
        sens_results = run_sensitivity_analyses(
            category_intervals, biser_by_chrom, chrom_sizes,
            eligible_by_chrom, sens_n, batch_size,
            np.random.default_rng(RANDOM_SEED + 1000),
        )

        # Write sensitivity results
        sens_path = TABLES_DIR / "purge_dups_biser_sensitivity_results.tsv"
        sens_headers = [
            "category", "analysis", "observed_overlap_bp", "observed_overlap_mb",
            "null_mean_overlap_bp", "null_sd_overlap_bp",
            "fold_enrichment", "empirical_p_enrichment", "n_permutations",
        ]
        with open(sens_path, "w") as f:
            f.write("\t".join(sens_headers) + "\n")
            for r in sens_results:
                f.write("\t".join(str(r.get(h, "")) for h in sens_headers) + "\n")
        print(f"\n  Wrote: {sens_path} ({len(sens_results)} sensitivity tests)")

    # Step 5: Chromosome-level analysis
    run_chromosome_analysis(category_intervals, biser_by_chrom, chrom_sizes)

    # Step 6: Runtime benchmark table
    bench_path = TABLES_DIR / "runtime_benchmark.tsv"
    with open(bench_path, "w") as f:
        f.write("run_mode\tn_permutations\tn_categories\tbatch_size\tn_workers\t"
                "elapsed_seconds\tpermutations_per_second\testimated_seconds_for_1000\t"
                "estimated_seconds_for_10000\n")
        f.write(f"{mode}\t{n_permutations}\t{len(CATEGORIES)}\t{batch_size}\t1\t"
                f"{bench['elapsed_seconds']:.1f}\t{bench['permutations_per_second']:.0f}\t"
                f"{bench['estimates'].get('pilot', 0):.0f}\t"
                f"{bench['estimates'].get('final', 0):.0f}\n")
    print(f"  Wrote: {bench_path}")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
