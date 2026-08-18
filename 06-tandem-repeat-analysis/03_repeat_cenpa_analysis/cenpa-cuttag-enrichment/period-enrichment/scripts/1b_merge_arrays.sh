#!/bin/bash
#SBATCH --job-name=trf_merge_arrays
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
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 1b_merge_arrays.sh
#
# Merge TRF intervals into repeat ARRAYS within each period bin, using
#   bedtools merge -d 0  (overlapping OR directly-bookended intervals merge;
#   NEVER across period bins). The merged arrays become the NEW unit of
#   observation for the period-enrichment analysis.
#
# This step is STATISTICAL PREP only — no signal is touched here. Signal is
# recalculated over merged arrays in 2_array_signal.sh (reusing the existing
# per-base bedGraph lookup tables — nothing is re-derived).
#
# New outputs (existing outputs are untouched):
#   data/merged/arrays/merged_arrays_all_bins.bed   chr,start,end,array_id,bin_id,n_intervals
#   data/merged/arrays/bin<N>_arrays.bed            per-bin: chr,start,end,array_id,n_intervals
#   data/merged/merged_array_statistics.csv         redundancy metrics per bin
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
# r-visualizations: this step needs Rscript (for the redundancy stats); bedtools
# is invoked by absolute path (BEDTOOLS_BIN), so no bulk-HiC-processing needed.
conda activate r-visualizations
set -euo pipefail

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
source "${SCRIPT_DIR}/config_period.sh"
init_period_dirs

BEDTOOLS_BIN="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/bulk-HiC-processing/bin/bedtools"
MERGED_DIR="${PERIOD_DATA_DIR}/merged"
ARRAY_DIR="${MERGED_DIR}/arrays"
mkdir -p "$ARRAY_DIR"

log "=== 1b: Merge TRF intervals into arrays (bedtools merge -d 0, per bin) ==="

# ── Step 1: merge per bin ─────────────────────────────────────────────────────
# For each bin: filter intervals to that bin (never across bins), sort by
# chr,start,end, then merge -d 0. `-c 4 -o count_distinct` counts how many
# distinct TRF interval_ids fell into each merged array (= intervals per array).
COMBINED_TMP="${ARRAY_DIR}/_merged_tmp.bed"
: > "${COMBINED_TMP}"
for b in $(seq 1 9); do
  awk -v b="$b" -F'\t' -v OFS='\t' '$5 == b {print $1, $2, $3, $4}' "${TRF_BED_CHR}" \
    | sort -k1,1 -k2,2n -k3,3n \
    | "${BEDTOOLS_BIN}" merge -d 0 -i stdin -c 4 -o count_distinct \
    | awk -v b="$b" -v OFS='\t' '{print $1, $2, $3, b, $4}' >> "${COMBINED_TMP}"
  log "  bin ${b}: $(awk -F'\t' -v b="$b" '$4 == b' "${COMBINED_TMP}" | wc -l) merged arrays"
done

# ── Step 2: assign globally-unique array_id; final column order ───────────────
# combined: chr, start, end, array_id, bin_id, n_intervals
COMBINED="${ARRAY_DIR}/merged_arrays_all_bins.bed"
awk -v OFS='\t' '{gid++; print $1, $2, $3, gid, $4, $5}' "${COMBINED_TMP}" > "${COMBINED}"
sort -k1,1 -k2,2n -k3,3n "${COMBINED}" -o "${COMBINED}"
rm -f "${COMBINED_TMP}"
log "  Total merged arrays across all bins: $(wc -l < "${COMBINED}")"

# ── Step 3: per-bin merged-array BEDs (shuffle inputs for the permutation) ────
for b in $(seq 1 9); do
  awk -v b="$b" -F'\t' -v OFS='\t' '$5 == b {print $1, $2, $3, $4, $6}' "${COMBINED}" \
    | sort -k1,1 -k2,2n -k3,3n > "${ARRAY_DIR}/bin${b}_arrays.bed"
  log "  bin${b}_arrays.bed: $(wc -l < "${ARRAY_DIR}/bin${b}_arrays.bed") arrays"
done

# ── Step 4: redundancy metrics per bin ────────────────────────────────────────
log "=== Computing redundancy metrics ==="
Rscript --no-save - "${COMBINED}" "${TRF_BED_CHR}" "${MERGED_DIR}/merged_array_statistics.csv" <<'REOF'
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(data.table))
merged <- fread(args[1], header = FALSE,
  col.names = c("chrom", "start", "end", "array_id", "bin_id", "n_intervals"))
trf    <- fread(args[2], header = FALSE,
  col.names = c("chrom", "start", "end", "interval_id", "bin_id",
                "period", "copies", "match_pct"))

bin_lab <- c("1"="1-10 bp","2"="11-50 bp","3"="51-192 bp","4"="193-195 bp","5"="196-347 bp",
             "6"="348-349 bp","7"="350-385 bp","8"="386-390 bp","9"="391+ bp")

int_stats <- trf[, .(n_intervals_total = .N, interval_bp = sum(end - start)), by = bin_id]
arr_stats <- merged[, .(
  n_arrays                   = .N,
  n_singleton_arrays         = sum(n_intervals == 1),
  n_multi_arrays             = sum(n_intervals >= 2),
  n_intervals_in_multi       = sum(n_intervals[n_intervals >= 2]),
  avg_intervals_per_array    = round(mean(n_intervals), 4),
  median_intervals_per_array = as.integer(median(n_intervals)),
  p75_intervals_per_array    = as.integer(quantile(n_intervals, .75)),
  max_intervals_per_array    = max(n_intervals),
  array_bp                   = sum(end - start)
), by = bin_id]

res <- merge(int_stats, arr_stats, by = "bin_id")
res[, bin_label := bin_lab[as.character(bin_id)]]
res[, frac_intervals_in_multi := round(n_intervals_in_multi / n_intervals_total, 4)]
res[, bp_overlap_frac := round((interval_bp - array_bp) / interval_bp, 4)]
setcolorder(res, c("bin_id", "bin_label", "n_intervals_total", "n_arrays",
                   "n_singleton_arrays", "n_multi_arrays", "n_intervals_in_multi",
                   "frac_intervals_in_multi", "avg_intervals_per_array",
                   "median_intervals_per_array", "p75_intervals_per_array",
                   "max_intervals_per_array", "interval_bp", "array_bp", "bp_overlap_frac"))
setorder(res, bin_id)
fwrite(res, args[3])
cat("\nRedundancy metrics (merged arrays, -d 0):\n")
print(res, digits = 4)
REOF
log "  Saved: ${MERGED_DIR}/merged_array_statistics.csv"

log "=== 1b complete ==="
