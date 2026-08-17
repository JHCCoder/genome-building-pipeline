#!/bin/bash
#SBATCH -J peaks2utr_subset_rerun
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
#SBATCH --no-requeue

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_PEAKS2UTR"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

# cwd matters for peaks2utr: CACHE_DIR = cwd/.cache, LOG_DIR = cwd/.log
cd "$PEAKS2UTR_DIR"

SUB="$SCRIPT_DIR/output/peaks2utr_changed_subset.gff3"   # from rebuild_subset.py
OUT="$SCRIPT_DIR/output/peaks2utr_subset_out.gff3"

peaks2utr "$SUB" merged.bam \
  --do-pseudo --keep-cache --extend-utr -p 8 --max-distance 1500 \
  -o "$OUT"

echo "peaks2utr subset finished; exit=$?"
