#!/bin/bash
#SBATCH -J peaks2utr_subset
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 02:00:00
#SBATCH --mem=32G
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
conda activate "$ENV_PEAKS2UTR"

# peaks2utr resolves .cache/.log relative to cwd, so run from its working
# directory (reuses cached stranded BAMs, MACS peaks and pileups).
cd "$PEAKS2UTR_WORK_DIR"

# Re-run peaks2utr on the changed-gene subset (built by 02_rebuild_changed_subset.py).
peaks2utr "$CHANGED_SUBSET_GFF" "$PEAKS2UTR_BAM" \
  --do-pseudo --keep-cache --extend-utr -p 8 --max-distance 1500 \
  -o "$PEAKS2UTR_SUBSET_OUT"

echo "peaks2utr subset finished; exit=$?"
