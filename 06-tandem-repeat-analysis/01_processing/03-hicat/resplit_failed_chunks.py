#!/usr/bin/env python3
"""
Re-split failed HiCAT chunk fastas into smaller sub-chunks to avoid OOM.

The HiCAT_HOR calculateED step builds an N×N edit distance matrix where
N = number of unique monomer blocks. Each failed chunk had:

  chr4:  ~181k blocks → N² ≈ 33B cells → ~262 GB matrix → OOM at 128G
  chr25: ~900k blocks → N² ≈ 730B cells → ~5.8 TB matrix → OOM at 400G

Target: ~25,000 blocks per sub-chunk → N² ≈ 625M cells → ~5 GB matrix.

Sub-chunks per failed chunk:
  chr4:  181k / 25k ≈ 8 sub-chunks  (density ~11.8 blocks/kb)
  chr25: 900k / 25k ≈ 36 sub-chunks (density ~73 blocks/kb)

Output structure:
  HiCAT_genome/<chr>/<chr>_partNN_subMM.fasta
  Header: ><chr>_partNN_subMM  (fixed later by run_HiCAT_subchunk.sh)

Provenance is preserved: chr4_part01 → chr4_part01_sub01..sub08
Merge step will concatenate sub-chunk outputs back to original chunk.
"""

import sys
from Bio import SeqIO
from pathlib import Path

BASE_DIR = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
                "/code/command-line-script/genome-annotation/HiCAT")
GENOME_DIR = BASE_DIR / "HiCAT_genome"
OVERLAP = 10_000  # 10 kb overlap each side (same as original split)

# ── Chunk definitions ──────────────────────────────────────────────
# Map: chromosome → list of (part_number, sub_chunks_to_create)
# Excludes chr4_part09 (job 11602606, still running with 500G)
# Target ~85-90k blocks per sub-chunk, anchored to known-good reference:
#   Job 7550536: full 10Mb chr4 centromere, ~86k blocks, 300G, 8 cores
#   Completed in 18h, used 159 GB actual memory.
#   chr4:  181k blocks / 90k ≈ 2 sub-chunks
#   chr25: 900k blocks / 90k ≈ 10 sub-chunks
FAILED_CHUNKS = {
    "chr4":  {"parts": ["01","02","03","04","05","06","07","08","10"],
              "sub_chunks": 2},
    "chr25": {"parts": ["01","02","03","04","05"],
              "sub_chunks": 10},
}


def split_chunk(chr_name: str, part: str, n_sub: int):
    """Split a single chunk fasta into n_sub sub-chunks with overlap."""
    chunk_fasta = GENOME_DIR / chr_name / f"{chr_name}_part{part}.fasta"

    if not chunk_fasta.exists():
        print(f"  SKIP: {chunk_fasta} not found")
        return

    rec = SeqIO.read(chunk_fasta, "fasta")
    total = len(rec.seq)
    sub_len = total // n_sub

    print(f"  {chr_name}_part{part}: {total:,} bp → {n_sub} sub-chunks "
          f"× ~{sub_len:,} bp (overlap {OVERLAP:,} bp)")

    for i in range(n_sub):
        raw_start = i * sub_len
        raw_end = (i + 1) * sub_len if i < n_sub - 1 else total

        start = max(0, raw_start - OVERLAP)
        end = min(total, raw_end + OVERLAP)

        sub_name = f"{chr_name}_part{part}_sub{i+1:02d}"
        sub_seq = rec.seq[start:end]
        sub_rec = rec.__class__(
            seq=sub_seq,
            id=sub_name,
            description=f"{chr_name}_part{part} chunk:{start+1}-{end}"
        )

        out_path = GENOME_DIR / chr_name / f"{sub_name}.fasta"
        SeqIO.write(sub_rec, out_path, "fasta")

    # Report summary
    sub_files = sorted((GENOME_DIR / chr_name).glob(f"{chr_name}_part{part}_sub*.fasta"))
    total_bp = sum(f.stat().st_size for f in sub_files)
    print(f"    → {len(sub_files)} sub-chunk fastas written "
          f"({total_bp:,} bp total)")


def main():
    dry_run = "--dry-run" in sys.argv

    for chr_name, cfg in FAILED_CHUNKS.items():
        print(f"\n{'='*60}")
        print(f"Processing {chr_name}: {len(cfg['parts'])} failed chunks, "
              f"{cfg['sub_chunks']} sub-chunks each "
              f"(~{len(cfg['parts']) * cfg['sub_chunks']} total sub-chunks)")
        print(f"{'='*60}")

        if dry_run:
            print("  [DRY RUN — no files written]")
            continue

        for part in cfg["parts"]:
            split_chunk(chr_name, part, cfg["sub_chunks"])

    print(f"\nDone. Sub-chunk fastas written to:")
    print(f"  {GENOME_DIR / 'chr4'}/chr4_part*_sub*.fasta")
    print(f"  {GENOME_DIR / 'chr25'}/chr25_part*_sub*.fasta")
    print(f"\nNext: run submit_subchunks.sh to submit the jobs.")


if __name__ == "__main__":
    main()
