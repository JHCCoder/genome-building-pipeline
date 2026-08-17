#!/bin/bash
#SBATCH -J biser_run
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 1-00:00:00
#SBATCH --mem=200G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Usage: sbatch 01_biser_run.sh <assembly> <species_tag>
#   <assembly>    genome fasta (full path, or a basename inside BISER_GENOME_DIR)
#   <species_tag> short label; becomes the output directory and BEDPE prefix

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

assembly="$1"
species="$2"
[ -n "$assembly" ] && [ -n "$species" ] || { echo "ERROR: usage: $0 <assembly> <species_tag>"; exit 1; }

# Resolve the assembly (full path, or basename within BISER_GENOME_DIR)
if [ ! -f "$assembly" ]; then
    assembly="$BISER_GENOME_DIR/$1"
fi

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_GENOME_ANNOTATION"

out_dir="$BISER_OUT_DIR/$species"
mkdir -p "$out_dir"
cd "$out_dir" || exit 1

# BISER expects a soft-masked genome; copy it into the output dir and index it
input=$(basename "$assembly")
cp "$assembly" .
samtools faidx "$input"

"$BISER_BIN" --threads "$BISER_THREADS" "$input" \
    --output "segdup_output_${species}" \
    --keep-contigs \
    --gc-heap "$BISER_GC_HEAP" \
    --max-edit-error "$BISER_MAX_EDIT_ERROR" \
    --max-error "$BISER_MAX_ERROR"
