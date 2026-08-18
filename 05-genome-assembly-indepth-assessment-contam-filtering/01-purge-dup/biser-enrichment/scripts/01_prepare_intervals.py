#!/usr/bin/env python3
"""
01_prepare_intervals.py — Parse, standardize, and validate purge_dups and BISER input files.

Reads the selected input files, produces standardized BED files, builds the BISER
nonredundant union, defines eligible genomic space, and generates summary tables
and validation reports.

Usage:
    python scripts/01_prepare_intervals.py
"""

import json
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import statistics

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj")
OUTPUT_ROOT = PROJECT_ROOT / "figure" / "segdup-purgeDup-overlap" / "purge-dups-biser-enrichment"
CONFIG_DIR = OUTPUT_ROOT / "config"
INTERMEDIATE_DIR = OUTPUT_ROOT / "intermediate"
TABLES_DIR = OUTPUT_ROOT / "tables"
LOGS_DIR = OUTPUT_ROOT / "logs"

# Load config
with open(CONFIG_DIR / "analysis_parameters.json") as f:
    CONFIG = json.load(f)

PURGE_DUPS_BED = Path(CONFIG["input_files"]["purge_dups_bed"])
BISER_BEDPE = Path(CONFIG["input_files"]["biser_bedpe"])
ASSEMBLY_FAI = Path(CONFIG["input_files"]["assembly_fai"])
CATEGORIES = CONFIG["analysis_parameters"]["categories"]
CATEGORY_ORDER = CATEGORIES
MITO_SEQ = CONFIG["analysis_parameters"]["mitochondrial_sequence"]


# ── Utility functions ──────────────────────────────────────────────────────

def run_cmd(cmd: List[str], desc: str = "") -> subprocess.CompletedProcess:
    """Run a subprocess command with error handling."""
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"ERROR running {desc}: {e}")
        print(f"STDERR: {e.stderr}")
        raise


def load_chromosome_sizes(fai_path: Path) -> Dict[str, int]:
    """Load chromosome sizes from FASTA index."""
    sizes: Dict[str, int] = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def validate_bed_coordinates(
    bed_path: Path,
    chrom_sizes: Dict[str, int],
    label: str = "",
) -> Tuple[List[str], List[str]]:
    """Validate BED coordinates. Returns (errors, warnings)."""
    errors: List[str] = []
    warnings: List[str] = []
    with open(bed_path) as f:
        for i, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                errors.append(f"{label} line {i}: fewer than 3 columns")
                continue
            chrom, start_str, end_str = parts[0], parts[1], parts[2]
            try:
                start, end = int(start_str), int(end_str)
            except ValueError:
                errors.append(f"{label} line {i}: non-integer coordinates")
                continue
            if start < 0:
                errors.append(f"{label} line {i}: negative start ({start})")
            if end <= start:
                errors.append(f"{label} line {i}: end ({end}) <= start ({start})")
            if chrom not in chrom_sizes:
                warnings.append(f"{label} line {i}: chromosome '{chrom}' not in assembly")
            elif end > chrom_sizes[chrom]:
                errors.append(
                    f"{label} line {i}: end ({end}) exceeds chromosome "
                    f"'{chrom}' length ({chrom_sizes[chrom]})"
                )
    return errors, warnings


def merge_intervals_python(intervals: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    """Merge overlapping or adjacent intervals within each chromosome using Python.

    Args:
        intervals: List of (chrom, start, end) tuples.

    Returns:
        Merged list of (chrom, start, end) tuples.
    """
    if not intervals:
        return []
    # Sort by chrom then start
    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    merged: List[Tuple[str, int, int]] = []
    for chrom, start, end in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            # Overlapping or adjacent: extend
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def write_bed(path: Path, intervals: List[Tuple[str, int, int]]) -> None:
    """Write BED3 intervals to file."""
    with open(path, "w") as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")


# ── Step 1: Load chromosome sizes ──────────────────────────────────────────

def step1_load_genome() -> Dict[str, int]:
    """Load chromosome sizes and write assembly.genome."""
    print("=" * 60)
    print("STEP 1: Loading chromosome sizes")
    chrom_sizes = load_chromosome_sizes(ASSEMBLY_FAI)
    total_bp = sum(chrom_sizes.values())
    print(f"  Assembly: {ASSEMBLY_FAI}")
    print(f"  Sequences: {len(chrom_sizes)}")
    print(f"  Total length: {total_bp:,} bp ({total_bp/1e6:.1f} Mb)")

    # Write assembly.genome if not exists
    genome_path = CONFIG_DIR / "assembly.genome"
    if not genome_path.exists():
        with open(genome_path, "w") as f:
            for chrom, size in chrom_sizes.items():
                f.write(f"{chrom}\t{size}\n")
        print(f"  Wrote: {genome_path}")

    # Print chromosome list
    chr_chroms = {k: v for k, v in chrom_sizes.items() if k.startswith("chr")}
    seq_chroms = {k: v for k, v in chrom_sizes.items() if k.startswith("seq")}
    other = {k: v for k, v in chrom_sizes.items() if not k.startswith(("chr", "seq"))}
    print(f"  Chromosome scaffolds (chr*): {len(chr_chroms)}")
    print(f"  Unplaced scaffolds (seq*): {len(seq_chroms)}")
    print(f"  Other sequences: {len(other)} ({list(other.keys())})")
    return chrom_sizes


# ── Step 2: Parse purge_dups categories ─────────────────────────────────────

def step2_parse_purge_dups(chrom_sizes: Dict[str, int]) -> Dict[str, List[Dict]]:
    """Parse purge_dups dups.bed into standardized category files.

    The dups.bed has 4 or 5 columns:
        chrom, start, end, category[, related_chrom]

    Returns dict mapping category -> list of interval dicts.
    """
    print("\n" + "=" * 60)
    print("STEP 2: Parsing purge_dups categories")

    category_intervals: Dict[str, List[Dict]] = {c: [] for c in CATEGORIES}
    all_intervals: List[Dict] = []
    seen = set()

    with open(PURGE_DUPS_BED) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                print(f"  WARNING: line {i} has {len(parts)} columns, skipping")
                continue

            chrom, start_str, end_str, category = parts[0], parts[1], parts[2], parts[3]
            start, end = int(start_str), int(end_str)
            related_chrom = parts[4] if len(parts) >= 5 else "."

            # Validate category
            if category not in CATEGORIES:
                print(f"  WARNING: line {i} has unknown category '{category}', skipping")
                continue

            # Validate coordinates
            if chrom not in chrom_sizes:
                print(f"  WARNING: line {i} chromosome '{chrom}' not in assembly, skipping")
                continue
            if start < 0 or end <= start or end > chrom_sizes[chrom]:
                print(f"  WARNING: line {i} invalid coordinates ({start}-{end}) for {chrom} (len={chrom_sizes[chrom]})")
                continue

            # Check for duplicates
            key = (chrom, start, end, category)
            if key in seen:
                print(f"  WARNING: line {i} duplicate interval: {chrom}:{start}-{end} {category}")
                continue
            seen.add(key)

            interval = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "category": category,
                "related_chrom": related_chrom,
                "length": end - start,
            }
            category_intervals[category].append(interval)
            all_intervals.append(interval)

    # Write standardized BED4
    bed4_path = INTERMEDIATE_DIR / "purge_dups.categories.bed"
    with open(bed4_path, "w") as f:
        for iv in sorted(all_intervals, key=lambda x: (x["chrom"], x["start"])):
            f.write(f"{iv['chrom']}\t{iv['start']}\t{iv['end']}\t{iv['category']}\n")
    print(f"  Wrote: {bed4_path} ({len(all_intervals)} intervals)")

    # Write per-category BED files
    for cat in CATEGORIES:
        cat_path = INTERMEDIATE_DIR / f"purge_dups.{cat}.bed"
        with open(cat_path, "w") as f:
            for iv in sorted(category_intervals[cat], key=lambda x: (x["chrom"], x["start"])):
                f.write(f"{iv['chrom']}\t{iv['start']}\t{iv['end']}\t{iv['category']}\n")
        print(f"  Wrote: {cat_path} ({len(category_intervals[cat])} intervals)")

    # Print summary
    print(f"\n  Category summary:")
    print(f"  {'Category':<12} {'N':>6} {'Total bp':>14} {'Total Mb':>12} {'Min':>10} {'Median':>10} {'Mean':>10} {'Max':>10}")
    print(f"  {'-'*12} {'-'*6} {'-'*14} {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for cat in CATEGORIES:
        ivs = category_intervals[cat]
        if not ivs:
            print(f"  {cat:<12} {'0':>6} {'0':>14} {'0.00':>12}")
            continue
        lengths = [iv["length"] for iv in ivs]
        total_bp = sum(lengths)
        print(
            f"  {cat:<12} {len(ivs):>6} {total_bp:>14,} {total_bp/1e6:>11.2f} "
            f"{min(lengths):>10,} {int(statistics.median(lengths)):>10,} "
            f"{statistics.mean(lengths):>10,.0f} {max(lengths):>10,}"
        )

    # Write category summary table
    summary_path = TABLES_DIR / "purge_dups_category_summary.tsv"
    with open(summary_path, "w") as f:
        f.write("category\tn_intervals\ttotal_bp\ttotal_mb\tmin_length\tmedian_length\tmean_length\tmax_length\n")
        for cat in CATEGORIES:
            ivs = category_intervals[cat]
            if not ivs:
                f.write(f"{cat}\t0\t0\t0.00\t0\t0\t0\t0\n")
                continue
            lengths = [iv["length"] for iv in ivs]
            f.write(
                f"{cat}\t{len(ivs)}\t{sum(lengths)}\t{sum(lengths)/1e6:.2f}\t"
                f"{min(lengths)}\t{int(statistics.median(lengths))}\t{statistics.mean(lengths):.0f}\t{max(lengths)}\n"
            )
    print(f"\n  Wrote: {summary_path}")

    # Chromosome distribution
    print(f"\n  Chromosome distribution by category:")
    for cat in CATEGORIES:
        chroms = set(iv["chrom"] for iv in category_intervals[cat])
        chr_count = sum(1 for c in chroms if c.startswith("chr"))
        seq_count = sum(1 for c in chroms if c.startswith("seq"))
        print(f"    {cat}: {len(chroms)} chromosomes ({chr_count} chr*, {seq_count} seq*)")

    # Check for cross-category overlaps using Python (avoid bedtools dependency)
    print(f"\n  Checking for within-category and cross-category overlaps...")
    # Sort intervals within each category by (chrom, start)
    for cat in CATEGORIES:
        category_intervals[cat].sort(key=lambda x: (x["chrom"], x["start"]))

    # Check within-category overlaps
    for cat in CATEGORIES:
        ivs = category_intervals[cat]
        for i in range(len(ivs) - 1):
            if ivs[i]["chrom"] == ivs[i+1]["chrom"] and ivs[i]["end"] > ivs[i+1]["start"]:
                print(f"  NOTE: Within-category overlap in {cat}: "
                      f"{ivs[i]['chrom']}:{ivs[i]['start']}-{ivs[i]['end']} vs "
                      f"{ivs[i+1]['chrom']}:{ivs[i+1]['start']}-{ivs[i+1]['end']}")

    # Check cross-category overlaps (simple pairwise)
    cross_overlaps = 0
    for i, cat1 in enumerate(CATEGORIES):
        for cat2 in CATEGORIES[i+1:]:
            for iv1 in category_intervals[cat1]:
                for iv2 in category_intervals[cat2]:
                    if iv1["chrom"] == iv2["chrom"] and iv1["start"] < iv2["end"] and iv2["start"] < iv1["end"]:
                        cross_overlaps += 1
    print(f"  Cross-category overlaps: {cross_overlaps}")

    return category_intervals


# ── Step 3: Parse BISER output ──────────────────────────────────────────────

def step3_parse_biser(chrom_sizes: Dict[str, int]) -> Tuple[Path, Path, int, int]:
    """Parse BISER BEDPE, extract both duplication arms, and build union.

    Returns (arms_bed_path, union_bed_path, n_pairs, n_arms).
    """
    print("\n" + "=" * 60)
    print("STEP 3: Parsing BISER output")

    arms_path = INTERMEDIATE_DIR / "biser.segmental_duplication_arms.bed"
    union_path = INTERMEDIATE_DIR / "biser.segmental_duplication_union.bed"

    n_pairs = 0
    n_arms = 0
    total_arm_bp = 0

    with open(BISER_BEDPE) as f_in, open(arms_path, "w") as f_out:
        for line in f_in:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            chrom1, start1, end1, chrom2, start2, end2 = parts[0:6]
            s1, e1, s2, e2 = int(start1), int(end1), int(start2), int(end2)

            # Validate
            if chrom1 in chrom_sizes and s1 >= 0 and e1 > s1 and e1 <= chrom_sizes[chrom1]:
                f_out.write(f"{chrom1}\t{s1}\t{e1}\n")
                n_arms += 1
                total_arm_bp += (e1 - s1)

            if chrom2 in chrom_sizes and s2 >= 0 and e2 > s2 and e2 <= chrom_sizes[chrom2]:
                f_out.write(f"{chrom2}\t{s2}\t{e2}\n")
                n_arms += 1
                total_arm_bp += (e2 - s2)

            n_pairs += 1

    print(f"  BISER pairs: {n_pairs:,}")
    print(f"  Extracted arms: {n_arms:,}")
    print(f"  Summed arm length: {total_arm_bp:,} bp ({total_arm_bp/1e6:.1f} Mb)")
    print(f"  Wrote: {arms_path}")

    # Sort and merge to create union using Python
    print(f"\n  Building nonredundant BISER union via Python merge...")
    arms_intervals: List[Tuple[str, int, int]] = []
    with open(arms_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                arms_intervals.append((parts[0], int(parts[1]), int(parts[2])))

    merged = merge_intervals_python(arms_intervals)
    write_bed(union_path, merged)

    # Count union stats
    n_union = 0
    union_bp = 0
    with open(union_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                n_union += 1
                union_bp += int(parts[2]) - int(parts[1])

    print(f"  Merged union intervals: {n_union:,}")
    print(f"  Nonredundant BISER length: {union_bp:,} bp ({union_bp/1e6:.1f} Mb)")
    print(f"  Wrote: {union_path}")

    # Genome-wide BISER fraction (excluding mito)
    total_genome = sum(v for k, v in chrom_sizes.items() if k != MITO_SEQ)
    genome_frac = union_bp / total_genome if total_genome > 0 else 0
    print(f"\n  Genome-wide BISER fraction (excl. mito): {genome_frac:.4f} ({genome_frac*100:.2f}%)")

    # Per-chromosome BISER summary
    chrom_biser: Dict[str, Dict] = defaultdict(lambda: {
        "n_intervals": 0, "biser_bp": 0, "chrom_length": 0
    })
    for chrom, size in chrom_sizes.items():
        chrom_biser[chrom]["chrom_length"] = size

    with open(union_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                if chrom in chrom_biser:
                    chrom_biser[chrom]["n_intervals"] += 1
                    chrom_biser[chrom]["biser_bp"] += (end - start)

    # Write chromosome summary
    chr_summary_path = TABLES_DIR / "biser_chromosome_summary.tsv"
    with open(chr_summary_path, "w") as f:
        f.write("chromosome\tchromosome_length\tn_biser_union_intervals\tbiser_union_bp\tbiser_union_mb\tbiser_fraction\n")
        # Sort: chr1-chr28, chrX, chrY, then seq* numerically
        def sort_key(item):
            name = item[0]
            if name.startswith("chr") and name[3:].isdigit():
                return (0, int(name[3:]))
            elif name == "chrX":
                return (0, 29)
            elif name == "chrY":
                return (0, 30)
            elif name.startswith("seq") and name[3:].isdigit():
                return (1, int(name[3:]))
            else:
                return (2, 0)

        for chrom, info in sorted(chrom_biser.items(), key=sort_key):
            frac = info["biser_bp"] / info["chrom_length"] if info["chrom_length"] > 0 else 0
            f.write(
                f"{chrom}\t{info['chrom_length']}\t{info['n_intervals']}\t"
                f"{info['biser_bp']}\t{info['biser_bp']/1e6:.2f}\t{frac:.6f}\n"
            )

    print(f"  Wrote: {chr_summary_path}")

    return arms_path, union_path, n_pairs, n_arms


# ── Step 4: Define eligible genomic space ──────────────────────────────────

def step4_eligible_genome(chrom_sizes: Dict[str, int]) -> Path:
    """Define eligible genomic space for random interval placement.

    Excludes mitochondrial sequence. All other sequences are eligible
    provided the interval fits within sequence boundaries.
    """
    print("\n" + "=" * 60)
    print("STEP 4: Defining eligible genomic space")

    eligible_path = INTERMEDIATE_DIR / "eligible_genome.bed"

    with open(eligible_path, "w") as f:
        for chrom, size in sorted(chrom_sizes.items()):
            if chrom == MITO_SEQ:
                continue
            f.write(f"{chrom}\t0\t{size}\n")

    total_eligible = sum(
        v for k, v in chrom_sizes.items() if k != MITO_SEQ
    )
    print(f"  Eligible sequences: {len(chrom_sizes) - 1}")
    print(f"  Eligible bases: {total_eligible:,} bp ({total_eligible/1e6:.1f} Mb)")
    print(f"  Excluded: {MITO_SEQ}")
    print(f"  Wrote: {eligible_path}")

    # Write eligibility summary
    summary_path = LOGS_DIR / "eligible_genome_summary.txt"
    with open(summary_path, "w") as f:
        f.write("Eligible Genome Summary\n")
        f.write("=======================\n\n")
        f.write(f"Total sequences in assembly: {len(chrom_sizes)}\n")
        f.write(f"Eligible sequences: {len(chrom_sizes) - 1}\n")
        f.write(f"Excluded: {MITO_SEQ}\n\n")
        f.write(f"Eligible bases: {total_eligible:,} bp ({total_eligible/1e6:.1f} Mb)\n\n")
        f.write("Per-chromosome eligible bases:\n")
        f.write(f"{'Chromosome':<15} {'Length (bp)':>15} {'Length (Mb)':>12}\n")
        f.write("-" * 42 + "\n")
        for chrom, size in sorted(chrom_sizes.items()):
            if chrom == MITO_SEQ:
                continue
            f.write(f"{chrom:<15} {size:>15,} {size/1e6:>11.2f}\n")

    print(f"  Wrote: {summary_path}")
    return eligible_path


# ── Step 5: Coordinate validation ──────────────────────────────────────────

def step5_validate_coordinates(
    chrom_sizes: Dict[str, int],
    category_intervals: Dict[str, List[Dict]],
    biser_union_path: Path,
) -> None:
    """Validate coordinate compatibility across all inputs."""
    print("\n" + "=" * 60)
    print("STEP 5: Validating coordinate compatibility")

    report_path = LOGS_DIR / "input_validation_report.txt"

    # Collect all names
    assembly_names = set(chrom_sizes.keys())
    purge_dups_names: Set[str] = set()
    for ivs in category_intervals.values():
        for iv in ivs:
            purge_dups_names.add(iv["chrom"])

    biser_names: Set[str] = set()
    with open(biser_union_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                biser_names.add(parts[0])

    shared = assembly_names & purge_dups_names & biser_names
    purge_only = purge_dups_names - assembly_names - biser_names
    biser_only = biser_names - assembly_names - purge_dups_names
    not_in_assembly = (purge_dups_names | biser_names) - assembly_names
    purge_not_in_biser = purge_dups_names - biser_names

    with open(report_path, "w") as f:
        f.write("Coordinate Compatibility Validation Report\n")
        f.write("==========================================\n\n")
        f.write(f"Assembly sequences: {len(assembly_names)}\n")
        f.write(f"purge_dups chromosomes: {len(purge_dups_names)}\n")
        f.write(f"BISER chromosomes: {len(biser_names)}\n")
        f.write(f"Shared by all three: {len(shared)}\n\n")

        if not_in_assembly:
            f.write(f"WARNING: Names not in assembly: {sorted(not_in_assembly)}\n\n")

        if purge_only:
            f.write(f"Names unique to purge_dups: {sorted(purge_only)}\n\n")
        if biser_only:
            f.write(f"Names unique to BISER: {sorted(biser_only)}\n\n")

        if purge_not_in_biser:
            f.write(f"purge_dups chromosomes absent from BISER:\n")
            for name in sorted(purge_not_in_biser):
                f.write(f"  {name}\n")
                # Count intervals on this chromosome
                for ivs in category_intervals.values():
                    for iv in ivs:
                        if iv["chrom"] == name:
                            f.write(f"    - {iv['category']}: {iv['chrom']}:{iv['start']}-{iv['end']} ({iv['length']} bp)\n")
            f.write("\n")

        # Report intervals on BISER-ineligible scaffolds
        f.write("Intervals on BISER-ineligible scaffolds:\n")
        total_ineligible_bp = 0
        for cat in CATEGORIES:
            for iv in category_intervals[cat]:
                if iv["chrom"] not in biser_names:
                    f.write(f"  {cat}: {iv['chrom']}:{iv['start']}-{iv['end']} ({iv['length']} bp)\n")
                    total_ineligible_bp += iv["length"]
        f.write(f"\nTotal category bases on BISER-ineligible scaffolds: {total_ineligible_bp:,} bp\n")

    print(f"  Shared names: {len(shared)}")
    print(f"  purge_dups chromosomes not in BISER: {len(purge_not_in_biser)}")
    if purge_not_in_biser:
        print(f"    Names: {sorted(purge_not_in_biser)}")
    print(f"  Wrote: {report_path}")


# ── Main ────────────────────────────────────────────────────────────────────

def main() -> None:
    """Run all preparation steps."""
    print("=" * 60)
    print("PURGE_DUPS vs BISER: Interval Preparation")
    print("=" * 60)
    print(f"purge_dups: {PURGE_DUPS_BED}")
    print(f"BISER:      {BISER_BEDPE}")
    print(f"Assembly:   {ASSEMBLY_FAI}")

    # Step 1: Load genome
    chrom_sizes = step1_load_genome()

    # Step 2: Parse purge_dups categories
    category_intervals = step2_parse_purge_dups(chrom_sizes)

    # Step 3: Parse BISER output
    arms_path, union_path, n_pairs, n_arms = step3_parse_biser(chrom_sizes)

    # Step 4: Define eligible genome
    eligible_path = step4_eligible_genome(chrom_sizes)

    # Step 5: Validate coordinates
    step5_validate_coordinates(chrom_sizes, category_intervals, union_path)

    print("\n" + "=" * 60)
    print("PREPARATION COMPLETE")
    print("=" * 60)
    print(f"Intermediate files in: {INTERMEDIATE_DIR}")
    print(f"Summary tables in:    {TABLES_DIR}")
    print(f"Validation logs in:   {LOGS_DIR}")


if __name__ == "__main__":
    main()
