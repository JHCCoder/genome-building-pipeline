#!/usr/bin/env python3
"""
Re-split chr25 sub-chunks into smaller a/b halves for rescue.

Each incomplete chr25 sub-chunk (~1.2 Mb, ~94k blocks) is split into
two sub-sub-chunks (~0.6 Mb, ~47k blocks each) to reduce ED phase
memory pressure and runtime.

  Original: 94k blocks → N² ≈ 8.8B comparisons → ~80 GB ED matrix
  Split:    47k blocks → N² ≈ 2.2B comparisons → ~18 GB ED matrix

Output:
  HiCAT_genome/chr25/chr25_partNN_subMMa.fasta
  HiCAT_genome/chr25/chr25_partNN_subMMb.fasta
  Header: >chr25_partNN_subMMa  (fixed later by run_HiCAT_rescue.sh)

Provenance: chr25_partNN → chr25_partNN_subMM → chr25_partNN_subMMa/b
Merge: concatenate a + b decomposition TSVs → original sub-chunk.

Usage:
  python resplit_chr25_rescue.py          # real run
  python resplit_chr25_rescue.py --dry-run
"""

import sys
from Bio import SeqIO
from pathlib import Path

BASE_DIR = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
                "/code/command-line-script/genome-annotation/HiCAT")
GENOME_DIR = BASE_DIR / "HiCAT_genome" / "chr25"
OVERLAP = 10_000  # 10 kb overlap each side

# Sub-chunks to split: all incomplete chr25 sub-chunks (36 total)
# These are sub-chunks that did NOT reach layer 12 with out/ dir.
# Determined 2026-07-29: part01_sub05-10, part02 all, part03 all, part04 all
SUBS_TO_SPLIT = [
    # part01: sub01-04 complete, sub05-10 need split
    ("01", ["05", "06", "07", "08", "09", "10"]),
    # part02: all 10 need split
    ("02", [f"{i:02d}" for i in range(1, 11)]),
    # part03: all 10 need split
    ("03", [f"{i:02d}" for i in range(1, 11)]),
    # part04: all 10 need split
    ("04", [f"{i:02d}" for i in range(1, 11)]),
]

N_SPLITS = 2  # split each sub-chunk into a and b halves


def split_subchunk(part: str, sub: str):
    """Split a single sub-chunk fasta into a/b halves with overlap."""
    sub_name = f"chr25_part{part}_sub{sub}"
    sub_fasta = GENOME_DIR / f"{sub_name}.fasta"

    if not sub_fasta.exists():
        print(f"  SKIP: {sub_fasta} not found")
        return

    rec = SeqIO.read(sub_fasta, "fasta")
    total = len(rec.seq)
    half = total // N_SPLITS

    print(f"  {sub_name}: {total:,} bp → {N_SPLITS} halves × ~{half:,} bp "
          f"(overlap {OVERLAP:,} bp)")

    suffixes = ["a", "b"]
    for i in range(N_SPLITS):
        raw_start = i * half
        raw_end = (i + 1) * half if i < N_SPLITS - 1 else total

        start = max(0, raw_start - OVERLAP)
        end = min(total, raw_end + OVERLAP)

        rescue_name = f"{sub_name}{suffixes[i]}"
        rescue_seq = rec.seq[start:end]
        rescue_rec = rec.__class__(
            seq=rescue_seq,
            id=rescue_name,
            description=f"{sub_name} chunk:{start+1}-{end}"
        )

        out_path = GENOME_DIR / f"{rescue_name}.fasta"
        SeqIO.write(rescue_rec, out_path, "fasta")
        print(f"    → {rescue_name}.fasta: {len(rescue_seq):,} bp "
              f"({start+1}-{end})")

    # Verify output
    a_fa = GENOME_DIR / f"{sub_name}a.fasta"
    b_fa = GENOME_DIR / f"{sub_name}b.fasta"
    a_len = len(SeqIO.read(a_fa, "fasta").seq) if a_fa.exists() else 0
    b_len = len(SeqIO.read(b_fa, "fasta").seq) if b_fa.exists() else 0
    print(f"    ✓ written: a={a_len:,} bp, b={b_len:,} bp")


def main():
    dry_run = "--dry-run" in sys.argv
    total = sum(len(subs) for _, subs in SUBS_TO_SPLIT)

    print(f"{'='*60}")
    print(f"chr25 rescue split: {total} sub-chunks → "
          f"{total * N_SPLITS} sub-sub-chunks "
          f"(~{total * N_SPLITS} new jobs)")
    print(f"{'='*60}")

    if dry_run:
        print("  [DRY RUN — no files written]")
        for part, subs in SUBS_TO_SPLIT:
            for sub in subs:
                sub_name = f"chr25_part{part}_sub{sub}"
                sub_fasta = GENOME_DIR / f"{sub_name}.fasta"
                if sub_fasta.exists():
                    rec = SeqIO.read(sub_fasta, "fasta")
                    half = len(rec.seq) // 2
                    print(f"  {sub_name}: {len(rec.seq):,} bp → "
                          f"{half:,} bp × 2")
        return

    count = 0
    for part, subs in SUBS_TO_SPLIT:
        print(f"\n-- part{part} ({len(subs)} sub-chunks) --")
        for sub in subs:
            split_subchunk(part, sub)
            count += 1

    print(f"\n{'='*60}")
    print(f"Done. {count} sub-chunks split into {count * N_SPLITS} "
          f"sub-sub-chunks.")
    print(f"Output: {GENOME_DIR}/chr25_part*_sub*[ab].fasta")
    print(f"\nNext: run submit_chr25_rescue.sh to submit the jobs.")


if __name__ == "__main__":
    main()
