#!/bin/bash
#SBATCH --job-name=trf_array_signal
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 2_array_signal.sh  (sbatch --array=0-3; one task per sample)
#
# Recalculate CENP-A CUT&Tag signal DIRECTLY over each merged repeat array
# (NOT by averaging the old interval-level values). Reuses the already-built
# run-length bedGraph lookup tables (data/permutation/bedgraph/<SAMPLE>.bedGraph)
# — the per-base genome signal is NOT re-derived.
#
# Output: data/merged/counts/trf_signal_array_<SAMPLE>.tsv
#   chr, start, end, array_id, bin_id, n_intervals, mean_coverage
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing
set -euo pipefail

source /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/config_period.sh
init_period_dirs

SAMPLES=("XG_150" "XG_151" "XG_152" "XG_153")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

ARRAY_DIR="${PERIOD_DATA_DIR}/merged/arrays"
COUNTS_DIR="${PERIOD_DATA_DIR}/merged/counts"
mkdir -p "$COUNTS_DIR"

COMBINED="${ARRAY_DIR}/merged_arrays_all_bins.bed"
BG="${PERIOD_DATA_DIR}/permutation/bedgraph/${SAMPLE}.bedGraph"
OUT="${COUNTS_DIR}/trf_signal_array_${SAMPLE}.tsv"

BEDTOOLS_BIN="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/bulk-HiC-processing/bin/bedtools"

log "=== 2_array_signal: ${SAMPLE} ==="
if [[ ! -f "$COMBINED" ]]; then
  log "ERROR: ${COMBINED} not found — run 1b_merge_arrays.sh first"
  exit 1
fi
if [[ ! -f "$BG" ]]; then
  log "ERROR: bedGraph not found: ${BG}"
  exit 1
fi

if [[ -f "$OUT" ]]; then
  log "  Exists: $(wc -l < "$OUT") rows — skip"
  exit 0
fi

log "  mapping ${SAMPLE} bedGraph over merged arrays..."
"${BEDTOOLS_BIN}" map -a "$COMBINED" -b "$BG" -c 4 -o mean -null 0 > "${OUT}.tmp"
mv "${OUT}.tmp" "$OUT"
log "  ${SAMPLE}: $(wc -l < "$OUT") arrays with mean coverage → ${OUT}"
log "=== 2_array_signal complete: ${SAMPLE} ==="
