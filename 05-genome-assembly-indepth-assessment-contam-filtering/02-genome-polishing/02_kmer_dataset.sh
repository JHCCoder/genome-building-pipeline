#!/bin/bash
#SBATCH -J polish_kmer
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 2-00:00:00
#SBATCH --mem=200G
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
conda activate "$ENV_GENOME_POLISHING"

mkdir -p "$POLISH_WORK_DIR"
cd "$POLISH_WORK_DIR"

# Build 21-mer and 31-mer datasets from short reads (remove -b 37 to count singletons)
yak count -t 8 -o male403_k21.yak -k 21 -b 37 "$POLISH_ILLUMINA_R1" "$POLISH_ILLUMINA_R2"
yak count -t 8 -o male403_k31.yak -k 31 -b 37 "$POLISH_ILLUMINA_R1" "$POLISH_ILLUMINA_R2"
