#!/bin/bash
#SBATCH -J busco
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 4:30:00
#SBATCH --mem=85G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p platinum
#SBATCH -q hcp-csd788
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_BUSCO"

prefix="haplotigPurged"
mkdir -p "$BUSCO_OUT_DIR"

# BUSCO in genome mode against the glires (rodent) lineage.
# The metrics (C/S/D/F/M, n) are in the resulting short_summary.*.txt file.
busco -i "$FINAL_ASSEMBLY" -m genome -l "$BUSCO_LINEAGE" -c 4 -f --out_path "$BUSCO_OUT_DIR" -o "$prefix"
