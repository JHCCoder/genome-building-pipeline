#!/bin/bash
#SBATCH -J purge_dup_step2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 3:00:00
#SBATCH --mem=12G
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

# Identify haplotigs/overlaps (dups.bed). Cutoffs/self-PAF are produced in step1.
"$PURGE_DUPS_BIN/purge_dups" -2 -T cutoffs -c PB.base.cov "$pri_asm.split.self.paf.gz" > dups.bed 2> purge_dups.log
