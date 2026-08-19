#!/bin/bash
#SBATCH -J peaks2utr_full
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 9:00:00
#SBATCH --mem=176G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p hotel
#SBATCH -q hotel
#SBATCH -A htl195
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

# Annotate 3' UTRs from RNA-seq coverage peaks (scMultiome GEX merged.bam).
# NOTE: peaks2utr generates ONLY 3' UTRs; 5' UTRs are carried through from the
# input GFF, never created by peaks2utr.
peaks2utr "$PEAKS2UTR_PREV_INPUT" "$PEAKS2UTR_BAM" \
  --do-pseudo --keep-cache --extend-utr -p 8 --max-distance 1500

echo "peaks2utr finished; exit=$?"
