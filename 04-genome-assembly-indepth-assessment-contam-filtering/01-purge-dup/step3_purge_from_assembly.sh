#!/bin/bash
#SBATCH -J purge_dup_step3
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --mem=128G
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
conda activate "$ENV_ASSESSMENT"

pri_asm="$CHR_ASSIGNED_ASSEMBLY"
cd "$PURGE_WORK_DIR"

# Remove the flagged dups from the assembly (produces hap.fa / purged.fa)
"$PURGE_DUPS_BIN/get_seqs" -e dups.bed "$pri_asm"
