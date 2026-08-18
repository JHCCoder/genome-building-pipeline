#!/bin/bash
#SBATCH -J salmon_merge
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 12:00:00
#SBATCH --mem=120G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc
conda activate "$ENV_TRANSCRIPTOME_MAPPING"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "$SCRIPT_DIR"

QUANT_DIRS=$(ls -d SRR*_salmon_quant)
echo "Found quantification directories:"
echo "$QUANT_DIRS"

salmon quantmerge \
  --quants $QUANT_DIRS \
  --names $(echo $QUANT_DIRS | sed 's/_salmon_quant//g') \
  --column tpm \
  --output merged_tpm_matrix.tsv
