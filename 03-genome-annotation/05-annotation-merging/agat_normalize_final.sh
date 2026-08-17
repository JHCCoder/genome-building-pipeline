#!/bin/bash
#SBATCH -J agat_normalize_final
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH --mem=32G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --no-requeue

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_AGAT"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "$SCRIPT_DIR"

# --- De-novo / degu-specific file paths (edit for your run) ---
INPUT="output/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.gff3"
OUTPUT="output/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.normalized.gff3"

# Normalize the assembled final annotation: rebuild gene(mRNA(CDS/exon/UTR/intron))
# Parent linkage (de-novo genes from braker_peak2utr.gff3 have orphaned mRNA/CDS/exon).
agat_convert_sp_gxf2gxf.pl -g "$INPUT" -o "$OUTPUT"

echo "AGAT exit=$?"
