#!/usr/bin/env python3
"""
Query per-base HiFi coverage for each gene locus (gene body +/- 5 kb flank).

For each chromosome, runs one `tabix` query covering the full chromosome, then
uses a sorted two-pointer merge to assign coverage intervals to regions (much
faster than per-region tabix calls).

Reads the region tables built by `01_build_regions.py` (regions.tsv and
busco_regions_5kb_flank.bed) and writes `region_coverage.tsv`.

Requires htslib's `tabix` (loaded via the `module load htslib/...` line below)
and a conda env with pandas/numpy.

    python 02_query_coverage.py [--max-genes N]   # --max-genes = quick test
"""

import argparse
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
# Project root — edit for your environment.
PROJ_ROOT = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"

# This script's directory (regions.tsv / busco beds are read from here;
# region_coverage.tsv is written here).
HERE = Path(__file__).resolve().parent

# Per-base HiFi coverage track (primary reads mapped to the final assembly).
COVERAGE_BED = (
    f"{PROJ_ROOT}/code/command-line-script/contig-coverage/"
    "hifiasm_041425_scaffolded_juiceBox_sorted_primaryReads_long_read_coverage_basepair_level.per-base.bed.gz"
)

# Chromosome lengths for the final assembly.
GENOME_LENGTH_FILE = (
    f"{PROJ_ROOT}/code/command-line-script/contig-coverage/"
    "hifiasm_041425_scaffolded_juiceBox_sorted_hardMasked_chrAssigned_genome_length.txt"
)

# htslib module setup (provides `tabix`). Adjust for your cluster/install.
HTSLIB_MODULES = (
    "module load shared cpu/0.21.2 gcc/13.3.0-c272c3y 2>/dev/null && "
    "module load htslib/1.17-f3pejh4 2>/dev/null"
)

# ============================================================
# Args
# ============================================================
parser = argparse.ArgumentParser()
parser.add_argument("--output", default=str(HERE / "region_coverage.tsv"))
parser.add_argument("--max-genes", type=int, default=None, help="Limit to first N genes for testing")
parser.add_argument("--chroms", type=str, default=None,
                    help="Comma-separated chromosomes to query (default: all in regions.tsv)")
args = parser.parse_args()

# ============================================================
# Load genome lengths
# ============================================================
print("Loading chromosome lengths...")
chr_lengths = {}
with open(GENOME_LENGTH_FILE) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            chr_lengths[parts[0]] = int(parts[1])

# ============================================================
# Load and prepare regions
# ============================================================
print("Loading regions...")
gene_meta = pd.read_csv(HERE / "regions.tsv", sep="\t")

regions_list = []
for _, row in gene_meta.iterrows():
    regions_list.append({
        "gene": row["gene"],
        "chrom": row["Chromosome"],
        "start": int(row["region_start"]),
        "end": int(row["region_end"]),
        "tier": row["tier"],
    })

regions = pd.DataFrame(regions_list)
print(f"  Gene regions: {len(regions)}")

busco_bed_path = HERE / "busco_regions_5kb_flank.bed"
if busco_bed_path.exists():
    busco = pd.read_csv(busco_bed_path, sep="\t", header=None,
                        names=["chrom", "start", "end", "gene", "tier", "overlaps_dup"])
    busco["start"] = busco["start"].astype(int)
    busco["end"] = busco["end"].astype(int)
    regions = pd.concat([regions, busco[["gene", "chrom", "start", "end", "tier"]]],
                        ignore_index=True)
    print(f"  With BUSCO: {len(regions)} total regions")

if args.max_genes:
    regions = regions.head(args.max_genes)
    print(f"  Limited to first {args.max_genes} regions (TEST MODE)")

if args.chroms:
    wanted = set(args.chroms.split(","))
    regions = regions[regions["chrom"].isin(wanted)].reset_index(drop=True)
    print(f"  Restricted to chromosomes {sorted(wanted)}: {len(regions)} regions")

regions = regions.sort_values(["chrom", "start"]).reset_index(drop=True)
regions["region_id"] = range(len(regions))

chr_to_regions = defaultdict(list)
for _, row in regions.iterrows():
    chr_to_regions[row["chrom"]].append(row)

CHROMOSOMES = sorted(
    chr_to_regions.keys(),
    key=lambda c: int(c.replace("chr", "")) if c.replace("chr", "").isdigit() else 999,
)
print(f"  Chromosomes: {CHROMOSOMES}")

# ============================================================
# For each chromosome, stream coverage and merge
# ============================================================
print(f"\nQuerying coverage for {len(regions)} regions across {len(CHROMOSOMES)} chromosomes...")


def weighted_median(intervals):
    """Weighted median from a list of (depth, length) tuples."""
    if not intervals:
        return np.nan
    sorted_int = sorted(intervals, key=lambda x: x[0])
    total = sum(l for _, l in sorted_int)
    half = total / 2
    cum = 0
    for depth, length in sorted_int:
        cum += length
        if cum >= half:
            return depth
    return sorted_int[-1][0]


results = []
start_time = time.time()
n_empty = 0

for chrom_idx, chrom in enumerate(CHROMOSOMES):
    chrom_regions = chr_to_regions[chrom]
    chrom_len = chr_lengths.get(chrom, 300_000_000)

    print(f"  [{chrom_idx + 1}/{len(CHROMOSOMES)}] {chrom}: {len(chrom_regions)} regions, "
          f"length={chrom_len / 1e6:.0f}Mb...", end=" ", flush=True)

    region_str = f"{chrom}:0-{chrom_len}"
    cmd = f"bash -c '{HTSLIB_MODULES} && tabix {COVERAGE_BED} {region_str}'"

    try:
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        print("TIMEOUT!")
        for r in chrom_regions:
            results.append({"gene": r["gene"], "tier": r["tier"], "chrom": chrom,
                            "mean_depth": np.nan, "median_depth": np.nan,
                            "region_length": r["end"] - r["start"]})
        continue

    cov_lines = proc.stdout.strip().split("\n")
    cov_intervals = []
    for line in cov_lines:
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 4:
            cov_intervals.append((int(parts[1]), int(parts[2]), int(parts[3])))

    t_chr_start = time.time()

    cov_idx = 0
    n_cov = len(cov_intervals)

    for region in chrom_regions:
        r_start = region["start"]
        r_end = region["end"]

        intervals = []
        total_bases = 0
        weighted_sum = 0.0
        weighted_sum_sq = 0.0
        zero_bases = 0
        min_depth = None
        max_depth = None

        while cov_idx < n_cov and cov_intervals[cov_idx][1] <= r_start:
            cov_idx += 1

        j = cov_idx
        while j < n_cov and cov_intervals[j][0] < r_end:
            c_start, c_end, depth = cov_intervals[j]
            overlap_start = max(r_start, c_start)
            overlap_end = min(r_end, c_end)
            overlap_len = overlap_end - overlap_start

            if overlap_len > 0:
                intervals.append((depth, overlap_len))
                total_bases += overlap_len
                weighted_sum += depth * overlap_len
                weighted_sum_sq += depth * depth * overlap_len
                if depth == 0:
                    zero_bases += overlap_len
                if min_depth is None or depth < min_depth:
                    min_depth = depth
                if max_depth is None or depth > max_depth:
                    max_depth = depth

            j += 1

        if total_bases == 0:
            n_empty += 1
            results.append({
                "gene": region["gene"], "tier": region["tier"], "chrom": chrom,
                "region_start": r_start, "region_end": r_end,
                "region_length": r_end - r_start,
                "bases_covered": 0, "mean_depth": np.nan,
                "median_depth": np.nan, "sum_depth": 0,
                "min_depth": np.nan, "max_depth": np.nan,
                "std_depth": np.nan, "pct_zero": np.nan,
            })
        else:
            mean_depth = weighted_sum / total_bases
            variance = (weighted_sum_sq / total_bases) - (mean_depth ** 2)
            std_depth = np.sqrt(max(0, variance))
            median_depth = weighted_median(intervals)
            pct_zero = zero_bases / total_bases * 100

            results.append({
                "gene": region["gene"], "tier": region["tier"], "chrom": chrom,
                "region_start": r_start, "region_end": r_end,
                "region_length": r_end - r_start,
                "bases_covered": total_bases, "mean_depth": mean_depth,
                "median_depth": median_depth, "sum_depth": weighted_sum,
                "min_depth": min_depth, "max_depth": max_depth,
                "std_depth": std_depth, "pct_zero": pct_zero,
            })

    t_elap = time.time() - t_chr_start
    print(f"{n_cov} cov intervals, {t_elap:.1f}s")

elapsed = time.time() - start_time
print(f"\nTotal: {len(results)} regions in {elapsed:.0f}s ({len(results) / elapsed:.1f}/s)")
print(f"Empty regions: {n_empty}")

# ============================================================
# Merge with gene metadata
# ============================================================
print("\nMerging with gene metadata...")
cov_df = pd.DataFrame(results)
cov_df = cov_df.merge(
    gene_meta[["gene", "base_gene", "num_copies", "gene_length", "Strand", "source", "is_paralog"]],
    on="gene", how="left",
)

cov_df.to_csv(args.output, sep="\t", index=False)
print(f"\nCoverage written to: {args.output}")
print(f"  Rows: {len(cov_df)}, Columns: {list(cov_df.columns)}")

print("\n" + "=" * 60)
print("COVERAGE SUMMARY BY TIER (mean depth)")
print("=" * 60)
for tier in ["T1", "T2", "T3", "T4", "BUSCO"]:
    subset = cov_df[cov_df["tier"] == tier]
    if len(subset) == 0:
        continue
    valid = subset["mean_depth"].dropna()
    print(f"  {tier:6s}: {len(subset):>6} regions  "
          f"median(mean)={valid.median():.2f}x  "
          f"mean(mean)={valid.mean():.2f}x  "
          f"NA={subset['mean_depth'].isna().sum()}")

print("\nDone!")
