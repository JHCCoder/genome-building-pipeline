#!/bin/bash
#SBATCH -J interproscan
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 06:00:00
#SBATCH --mem=30G
#SBATCH --array=66    # straggler rerun (single chunk); for a full run use --array=1-100
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

# Strip '*' (stop codons) from the chunk, then run InterProScan
sed 's/\*//g' 041425-braker-split-files-100files/chunk_${SLURM_ARRAY_TASK_ID}.fa \
  > 041425-braker-interproscan-result/chunk_${SLURM_ARRAY_TASK_ID}.clean.aa

"$INTERPROSCAN_BIN" \
  -i 041425-braker-interproscan-result/chunk_${SLURM_ARRAY_TASK_ID}.clean.aa \
  -f tsv -dp \
  -o 041425-braker-interproscan-result/interproscan_result_${SLURM_ARRAY_TASK_ID} \
  -cpu 4
