#!/usr/bin/env python3
"""
Prepare HiCAT decomposition data for karyoploteR visualization.

Parses HiCAT final_decomposition.tsv files, bins monomers by genomic
position, and outputs BED files for plotting with karyoploteR.

HiCAT output columns (tab-separated, no header):
  1. chromosome
  2. monomer_id (format: chr:start:end_type, e.g. chr4:82735000:82736000_1)
  3. array_start_idx  - monomer start index within HOR decomposition
  4. array_end_idx    - monomer end index within HOR decomposition
  5. identity_pct     - percent identity to monomer consensus
  6-11. extra metrics (mostly None / -1.00)
  12. strand

Output files:
  - hicat_{chr}_density.bed  : binned monomer density (count per bin)
  - hicat_{chr}_types.bed    : dominant monomer type per bin
"""

import sys
import os
from collections import defaultdict, Counter

# ── Configuration ────────────────────────────────────────────────────────────

HICAT_DIR = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/genome-annotation/HiCAT/HiCAT_genome"
OUT_DIR = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"

# Chromosome sizes from FASTA index
CHR_SIZES = {
    "chr4":  151562733,
    "chr25": 60231632,
}

# Chromosomes to process
CHROMOSOMES = ["chr4", "chr25"]

# Bin size for density track
BIN_SIZE = 200000  # 200 kb

# Gap threshold for identifying separate HOR arrays (bp)
# Monomers separated by more than this are considered different arrays
HOR_GAP_THRESHOLD = 50000  # 50 kb


def parse_monomer_id(monomer_id):
    """
    Parse a HiCAT monomer ID like 'chr4:82735000:82736000_1'
    Returns (genomic_start, genomic_end, monomer_type)
    """
    # Format: chr4:82735000:82736000_1
    parts = monomer_id.split(":")
    # parts[0] = chr, parts[1] = start, parts[2] = "end_type"
    gen_start = int(parts[1])
    rest = parts[2]
    # Split on last underscore to separate end from type
    # Type may contain quotes like "1'" so find last underscore
    underscore_idx = rest.rindex("_")
    gen_end = int(rest[:underscore_idx])
    monomer_type = rest[underscore_idx + 1:]
    return gen_start, gen_end, monomer_type


def process_chromosome(chrom):
    """
    Parse HiCAT decomposition for one chromosome.
    Returns:
      bins: dict mapping bin_start -> Counter of monomer types
      monomers: list of (gen_start, gen_end, monomer_type) sorted by position
    """
    tsv_path = os.path.join(HICAT_DIR, chrom, "final_decomposition.tsv")
    if not os.path.exists(tsv_path):
        print(f"ERROR: File not found: {tsv_path}", file=sys.stderr)
        return None, None

    chr_size = CHR_SIZES[chrom]
    n_bins = (chr_size // BIN_SIZE) + 1

    # Initialize bins
    bins = defaultdict(Counter)  # bin_start -> Counter of monomer types
    total_monomers = 0

    print(f"  Reading {tsv_path} ...", file=sys.stderr)

    with open(tsv_path) as f:
        for line_i, line in enumerate(f):
            if line_i % 1000000 == 0 and line_i > 0:
                print(f"    {line_i/1e6:.0f}M lines processed, {total_monomers} monomers binned ...",
                      file=sys.stderr)

            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            monomer_id = parts[1]
            try:
                gen_start, gen_end, monomer_type = parse_monomer_id(monomer_id)
            except (ValueError, IndexError) as e:
                if line_i < 5:
                    print(f"    WARNING: could not parse monomer_id '{monomer_id}': {e}",
                          file=sys.stderr)
                continue

            # Determine bin
            bin_idx = gen_start // BIN_SIZE
            bin_start = bin_idx * BIN_SIZE

            bins[bin_start][monomer_type] += 1
            total_monomers += 1

    print(f"    Done. {total_monomers:,} monomers binned into {len(bins)} bins.",
          file=sys.stderr)
    return bins, total_monomers


def write_density_bed(chrom, bins, out_dir):
    """Write monomer density BED file."""
    out_path = os.path.join(out_dir, f"hicat_{chrom}_density.bed")
    with open(out_path, "w") as f:
        f.write(f"# HiCAT monomer density for {chrom}\n")
        f.write(f"# Bin size: {BIN_SIZE} bp\n")
        f.write(f"# Columns: chrom, start, end, monomer_count, density_per_kb\n")
        for bin_start in sorted(bins.keys()):
            total = sum(bins[bin_start].values())
            bin_end = min(bin_start + BIN_SIZE, CHR_SIZES[chrom])
            density_per_kb = total / ((bin_end - bin_start) / 1000.0) if bin_end > bin_start else 0
            f.write(f"{chrom}\t{bin_start}\t{bin_end}\t{total}\t{density_per_kb:.2f}\n")
    print(f"  Wrote: {out_path}", file=sys.stderr)
    return out_path


def write_types_bed(chrom, bins, out_dir):
    """Write dominant monomer type per bin BED file."""
    out_path = os.path.join(out_dir, f"hicat_{chrom}_types.bed")
    with open(out_path, "w") as f:
        f.write(f"# HiCAT dominant monomer type per bin for {chrom}\n")
        f.write(f"# Bin size: {BIN_SIZE} bp\n")
        f.write(f"# Columns: chrom, start, end, dominant_type, type_count, total_monomers, purity\n")
        for bin_start in sorted(bins.keys()):
            counter = bins[bin_start]
            total = sum(counter.values())
            if total == 0:
                continue
            dominant_type, dominant_count = counter.most_common(1)[0]
            purity = dominant_count / total if total > 0 else 0
            bin_end = min(bin_start + BIN_SIZE, CHR_SIZES[chrom])
            f.write(f"{chrom}\t{bin_start}\t{bin_end}\t{dominant_type}\t{dominant_count}\t{total}\t{purity:.3f}\n")
    print(f"  Wrote: {out_path}", file=sys.stderr)
    return out_path


def write_hor_arrays_bed(chrom, out_dir):
    """
    Identify contiguous HOR arrays by sorting monomers by genomic position
    and merging those within HOR_GAP_THRESHOLD. Output as BED file.

    A HOR array is a contiguous stretch of monomers.
    """
    tsv_path = os.path.join(HICAT_DIR, chrom, "final_decomposition.tsv")

    print(f"  Identifying HOR arrays for {chrom} ...", file=sys.stderr)

    # Collect all monomer positions
    positions = []  # list of (gen_start, gen_end, monomer_type)
    with open(tsv_path) as f:
        for line_i, line in enumerate(f):
            if line_i % 1000000 == 0 and line_i > 0:
                print(f"    {line_i/1e6:.0f}M lines ...", file=sys.stderr)
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            try:
                gen_start, gen_end, monomer_type = parse_monomer_id(parts[1])
                positions.append((gen_start, gen_end, monomer_type))
            except (ValueError, IndexError):
                continue

    print(f"    Sorting {len(positions):,} monomers by genomic position ...", file=sys.stderr)
    positions.sort(key=lambda x: x[0])

    # Merge into arrays
    arrays = []
    if positions:
        arr_start = positions[0][0]
        arr_end = positions[0][1]
        arr_types = Counter([positions[0][2]])
        arr_n = 1

        for gen_start, gen_end, mtype in positions[1:]:
            if gen_start - arr_end > HOR_GAP_THRESHOLD:
                # Finish current array
                arrays.append({
                    "start": arr_start,
                    "end": arr_end,
                    "n_monomers": arr_n,
                    "dominant_type": arr_types.most_common(1)[0][0],
                    "n_types": len(arr_types),
                })
                arr_start = gen_start
                arr_end = gen_end
                arr_types = Counter([mtype])
                arr_n = 1
            else:
                arr_end = max(arr_end, gen_end)
                arr_types[mtype] += 1
                arr_n += 1

        # Don't forget last array
        arrays.append({
            "start": arr_start,
            "end": arr_end,
            "n_monomers": arr_n,
            "dominant_type": arr_types.most_common(1)[0][0],
            "n_types": len(arr_types),
        })

    # Filter: only keep arrays with >= 2 monomers (singletons aren't HORs)
    arrays = [a for a in arrays if a["n_monomers"] >= 2]

    print(f"    Found {len(arrays):,} HOR arrays (>=2 monomers, gap <= {HOR_GAP_THRESHOLD} bp)",
          file=sys.stderr)

    out_path = os.path.join(out_dir, f"hicat_{chrom}_hor_arrays.bed")
    with open(out_path, "w") as f:
        f.write(f"# HiCAT HOR arrays for {chrom}\n")
        f.write(f"# Gap threshold: {HOR_GAP_THRESHOLD} bp\n")
        f.write(f"# Columns: chrom, start, end, n_monomers, dominant_type, n_types\n")
        for a in arrays:
            f.write(f"{chrom}\t{a['start']}\t{a['end']}\t"
                    f"{a['n_monomers']}\t{a['dominant_type']}\t{a['n_types']}\n")
    print(f"  Wrote: {out_path}", file=sys.stderr)

    # Print summary stats
    sizes = [a["end"] - a["start"] + 1 for a in arrays]
    if sizes:
        print(f"    Array size range: {min(sizes):,} - {max(sizes):,} bp", file=sys.stderr)
        print(f"    Array size median: {sorted(sizes)[len(sizes)//2]:,} bp", file=sys.stderr)
        total_bp = sum(sizes)
        print(f"    Total array coverage: {total_bp:,} bp ({100*total_bp/CHR_SIZES[chrom]:.1f}% of chromosome)",
              file=sys.stderr)

    return out_path


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    for chrom in CHROMOSOMES:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"Processing {chrom} ({CHR_SIZES[chrom]:,} bp)", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

        # Parse and bin monomers
        bins, total = process_chromosome(chrom)
        if bins is None:
            continue

        # Write density BED
        write_density_bed(chrom, bins, OUT_DIR)

        # Write dominant type BED
        write_types_bed(chrom, bins, OUT_DIR)

        # Write HOR arrays BED
        write_hor_arrays_bed(chrom, OUT_DIR)

    print(f"\nDone. Output files in: {OUT_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
