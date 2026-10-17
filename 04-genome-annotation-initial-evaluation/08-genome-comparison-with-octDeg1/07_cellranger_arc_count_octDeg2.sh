#!/bin/bash
#SBATCH -J 101226_cellrangerArc_count_octDeg2_liftoffOnly
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 36:00:00
#SBATCH --mem=128G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue

# Re-run Cell Ranger ARC count for ONE scMultiome sample on OctDeg2.0 with the Liftoff-only
# reference ($FIG3D_OCTDEG2_CR_REF, built by 06_cellranger_arc_mkref_liftoff.sh), so Figure 3D
# compares OctDeg1.0 and OctDeg2.0 on matched gene models.
#
# The libraries per sample come from data/fig3d_multiome_libraries.tsv and are the same
# FASTQ sets as the original OctDeg1 runs (column octdeg1_run). This corrects the July 2025
# OctDeg2 runs, where the run named 181_dHIP used 6997_dHIP's libraries (WB_173/WB_174).
#
# Usage (one job per sample):
#   sbatch --mail-user=<you> 07_cellranger_arc_count_octDeg2.sh 181_dHIP

set -euo pipefail

SAMPLE="${1:?usage: 07_cellranger_arc_count_octDeg2.sh <sample from data/fig3d_multiome_libraries.tsv>}"

# Under sbatch, BASH_SOURCE is Slurm's spool copy of this script; prefer SLURM_SUBMIT_DIR.
if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    _script_dir="$SLURM_SUBMIT_DIR"
else
    _script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
fi

# Load shared configuration (config.sh at the repo root)
_repo_root="$_script_dir"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

LIB_TSV="$_script_dir/data/fig3d_multiome_libraries.tsv"
row="$(awk -F'\t' -v s="$SAMPLE" 'NR>1 && $1==s' "$LIB_TSV")"
[[ -n "$row" ]] || { echo "ERROR: sample $SAMPLE not in $LIB_TSV" >&2; exit 1; }
IFS=$'\t' read -r _ ATAC GEX _ <<< "$row"

[[ -f "$FIG3D_OCTDEG2_CR_REF/reference.json" ]] || { echo "ERROR: run 06_cellranger_arc_mkref_liftoff.sh first" >&2; exit 1; }

RUN_ID="${SAMPLE}_OctDeg2_liftoffOnly"
COUNT_DIR="$FIG3D_CR_DIR/count"
mkdir -p "$COUNT_DIR"
cd "$COUNT_DIR"                 # cellranger-arc writes ./$RUN_ID
if [[ -e "$RUN_ID" ]]; then
  echo "ERROR: $COUNT_DIR/$RUN_ID already exists; move it aside before re-running" >&2
  exit 1
fi

# Library CSV; every listed library must have FASTQs in the readable Lustre FASTQ directory.
LIB_CSV="$COUNT_DIR/${RUN_ID}.libraries.csv"
{
  echo "fastqs,sample,library_type"
  for lib in ${ATAC//,/ }; do echo "$FIG3D_MULTIOME_FASTQ_DIR,$lib,Chromatin Accessibility"; done
  for lib in ${GEX//,/ };  do echo "$FIG3D_MULTIOME_FASTQ_DIR,$lib,Gene Expression"; done
} > "$LIB_CSV"
for lib in ${ATAC//,/ } ${GEX//,/ }; do
  ls "$FIG3D_MULTIOME_FASTQ_DIR/${lib}"_S*_R*_001.fastq.gz >/dev/null 2>&1 \
    || { echo "ERROR: no FASTQs for $lib in $FIG3D_MULTIOME_FASTQ_DIR" >&2; exit 1; }
done
cat "$LIB_CSV"

echo "cellranger-arc: $("$CELLRANGER_ARC_BIN" --version 2>&1 | head -1)"
"$CELLRANGER_ARC_BIN" count --id="$RUN_ID" \
  --reference="$FIG3D_OCTDEG2_CR_REF" \
  --libraries="$LIB_CSV" \
  --localcores=16 --localmem=120

echo "count finished: $(date)"
