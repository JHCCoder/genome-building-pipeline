#!/usr/bin/env python3
"""Split chr4.fasta into N equal-length chunks with a small overlap for boundary-spanning monomers."""
import sys
from Bio import SeqIO
from pathlib import Path

IN_FASTA = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/genome-annotation/HiCAT/HiCAT_genome/chr4/chr4.fasta")
OUT_DIR  = IN_FASTA.parent
N_CHUNKS = 10
OVERLAP  = 10_000  # 10 kb overlap on each side to avoid missing boundary monomers

def main():
    rec = SeqIO.read(IN_FASTA, "fasta")
    total = len(rec.seq)
    chunk_len = total // N_CHUNKS
    print(f"chr4 total length: {total:,} bp")
    print(f"Chunks: {N_CHUNKS} × ~{chunk_len:,} bp with {OVERLAP:,} bp overlap\n")

    for i in range(N_CHUNKS):
        raw_start = i * chunk_len          # 0-based
        raw_end   = (i + 1) * chunk_len if i < N_CHUNKS - 1 else total

        # Apply overlap
        start = max(0, raw_start - OVERLAP)
        end   = min(total, raw_end + OVERLAP)

        part_name = f"chr4_part{i+1:02d}"
        chunk_seq = rec.seq[start:end]
        chunk_rec = rec.__class__(seq=chunk_seq, id=part_name, description=f"chr4:{start+1}-{end}")

        out_path = OUT_DIR / f"{part_name}.fasta"
        SeqIO.write(chunk_rec, out_path, "fasta")
        print(f"  {part_name}: {start+1:,}-{end:,} ({end-start:,} bp) → {out_path}")

if __name__ == "__main__":
    main()
