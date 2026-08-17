#!/bin/bash
#SBATCH -J merqury
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 4:00:00
#SBATCH --mem=48G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_MERQURY"

prefix="hifiasm_0414_primary_contamR_haplotigR"
mkdir -p "$MERQURY_OUT_DIR/$prefix"
cd "$MERQURY_OUT_DIR/$prefix"

# Link (or build via 01_meryl_db.sh) the Illumina WGS k-mer database
ln -sf "$MERQURY_OUT_DIR/$MERYL_DB" .

# Assess assembly completeness/consensus vs. the read k-mer db.
# Usage: merqury.sh <read.meryl> <assembly.fasta> <output_prefix>
merqury.sh "$MERYL_DB" "$FINAL_ASSEMBLY" "$prefix"
