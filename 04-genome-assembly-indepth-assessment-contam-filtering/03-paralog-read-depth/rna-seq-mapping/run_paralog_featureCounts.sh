#!/bin/bash
#SBATCH -J paralog_featureCounts
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 2-00:00:00
#SBATCH --mem=150G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue

# Quantify RNA-seq read support for all paralog families
# (2,252 paralogs + 1,095 parents = 3,347 gene features) across 29 tissue BAMs.

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc
conda activate "$ENV_TRANSCRIPTOME_MAPPING"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "$SCRIPT_DIR"

ANNOTATION="paralog_families.gff"   # from build_paralog_families.py
OUTPUT="paralog_families_counts.tsv"

echo "Gene features: $(awk -F'\t' '$3=="gene"' "$ANNOTATION" | wc -l)"

BAM_FILES=( *_Aligned.sortedByCoord.out.bam )
echo "BAM files found: ${#BAM_FILES[@]}"

featureCounts \
    -a "$ANNOTATION" \
    -o "$OUTPUT" \
    -g ID \
    -t gene \
    -p \
    -T 16 \
    --primary \
    "${BAM_FILES[@]}"

echo "featureCounts exit code: $?"
echo "Finished at: $(date)"
