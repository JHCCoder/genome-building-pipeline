#!/usr/bin/env python3
"""Split a chromosome fasta into N equal-length chunks with overlap for boundary-spanning monomers.

Usage: python split_chromosome.py <chr> <N_chunks> [overlap_bp]

Outputs chunks to HiCAT_genome/<chr>/<chr>_partNN.fasta
"""
import sys
from Bio import SeqIO
from pathlib import Path

BASE_DIR = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/genome-annotation/HiCAT")
ASSEMBLY = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta")
SAMTOOLS = "/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-samtools/bin/samtools"

def main():
    chr_name = sys.argv[1]
    n_chunks = int(sys.argv[2])
    overlap  = int(sys.argv[3]) if len(sys.argv) > 3 else 10_000

    out_dir = BASE_DIR / "HiCAT_genome" / chr_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # Extract chromosome from assembly
    chr_fasta = out_dir / f"{chr_name}.fasta"
    if not chr_fasta.exists() or chr_fasta.stat().st_size == 0:
        import subprocess
        print(f"Extracting {chr_name} from assembly...")
        result = subprocess.run(
            [SAMTOOLS, "faidx", str(ASSEMBLY), chr_name],
            capture_output=True, text=True
        )
        with open(chr_fasta, 'w') as f:
            f.write(result.stdout)
        print(f"  Extracted to {chr_fasta}")

    rec = SeqIO.read(chr_fasta, "fasta")
    total = len(rec.seq)
    chunk_len = total // n_chunks
    print(f"{chr_name} total length: {total:,} bp")
    print(f"Chunks: {n_chunks} × ~{chunk_len:,} bp with {overlap:,} bp overlap\n")

    for i in range(n_chunks):
        raw_start = i * chunk_len
        raw_end   = (i + 1) * chunk_len if i < n_chunks - 1 else total

        start = max(0, raw_start - overlap)
        end   = min(total, raw_end + overlap)

        part_name = f"{chr_name}_part{i+1:02d}"
        chunk_seq = rec.seq[start:end]
        chunk_rec = rec.__class__(seq=chunk_seq, id=part_name, description=f"{chr_name}:{start+1}-{end}")

        out_path = out_dir / f"{part_name}.fasta"
        SeqIO.write(chunk_rec, out_path, "fasta")
        print(f"  {part_name}: {start+1:,}-{end:,} ({end-start:,} bp)")

if __name__ == "__main__":
    main()
