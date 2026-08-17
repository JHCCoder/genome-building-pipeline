#!/bin/bash
#SBATCH -J discoverY
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 12:00:00
#SBATCH --mem=360G
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
conda activate "$ENV_DISCOVERY"

# Optional bloom-filter mode (faster for large datasets):
#python discoverY.py --female_kmers_set --female_bloom --kmer_size 25 --mode female+male --female_bloom_capacity 3100000000

python discoverY.py --female_kmers_set --kmer_size 25 --mode female+male
