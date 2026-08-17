#!/bin/bash
#SBATCH --job-name=trf_bedgraph
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 2e_prepare_bedgraphs.sh
# Convert each sample's per-base CUT&Tag coverage (8 GB gz) to a run-length
# encoded bedGraph ONCE, into a shared location. All per-bin permutation jobs
# (2d) read these, so the conversion is never repeated per bin.
#
# Array: 1 task per sample (XG_150/151/152/153).
# Output: data/permutation/bedgraph/<SAMPLE>.bedGraph
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing
set -euo pipefail

source /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/config_period.sh
init_period_dirs

SAMPLES=("XG_150" "XG_151" "XG_152" "XG_153")
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
COV_FILE="${COVERAGE_DIR}/${SAMPLE}_perbase.txt.gz"
BG_DIR="${PERIOD_DATA_DIR}/permutation/bedgraph"
mkdir -p "$BG_DIR"
BG_FILE="${BG_DIR}/${SAMPLE}.bedGraph"

log "=== BedGraph prep for ${SAMPLE} ==="

if [[ -f "${BG_FILE}" ]]; then
  N_BAD=$(awk -F'\t' 'NF!=4 {n++} END {print n+0}' "${BG_FILE}")
  if [[ "$N_BAD" -eq 0 ]]; then
    log "  Exists (valid): $(wc -l < ${BG_FILE}) blocks — skip"
    exit 0
  fi
  log "  Existing bedGraph corrupt (${N_BAD} bad lines) — regenerating"
  rm -f "${BG_FILE}"
fi

log "  Converting per-base coverage to bedGraph..."
zcat "${COV_FILE}" 2>/dev/null | awk '
BEGIN { FS="\t"; OFS="\t" }
NR==1 {
    chr=$1; pos=$2; val=$3
    run_start=pos-1; run_end=pos; run_chr=chr; run_val=val
    next
}
{
    if ($1 != run_chr || $3 != run_val) {
        print run_chr, run_start, run_end, run_val
        run_chr=$1; run_start=$2-1; run_end=$2; run_val=$3
    } else {
        run_end=$2
    }
}
END {
    if (NR > 0) print run_chr, run_start, run_end, run_val
}' > "${BG_FILE}"

# Validate + truncate any trailing corrupt lines from a truncated decompression
N_BAD=$(awk -F'\t' 'NF!=4 {n++} END {print n+0}' "${BG_FILE}")
if [[ "$N_BAD" -gt 0 ]]; then
  FIRST_BAD=$(awk -F'\t' 'NF!=4 {print NR; exit}' "${BG_FILE}")
  log "  WARNING: ${N_BAD} corrupt line(s) — truncating at $((FIRST_BAD - 1))"
  head -n $((FIRST_BAD - 1)) "${BG_FILE}" > "${BG_FILE}.clean"
  mv "${BG_FILE}.clean" "${BG_FILE}"
fi

N_BLOCKS=$(wc -l < "${BG_FILE}")
N_BAD=$(awk -F'\t' 'NF!=4 {n++} END {print n+0}' "${BG_FILE}")
log "  ${SAMPLE}: ${N_BLOCKS} bedGraph blocks, ${N_BAD} bad lines"
log "  Saved: ${BG_FILE}"
log "=== BedGraph prep complete: ${SAMPLE} ==="
