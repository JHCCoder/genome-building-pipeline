#!/bin/bash
#SBATCH --job-name=trf_combine_analyze
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 2f_combine_and_analyze.sh
# Runs AFTER all per-bin permutation jobs (2d, array 1-9) complete.
#
#   1. Combine per-bin ctrl nulls (bins 1-5) -> data/counts/trf_ctrl_bg_signal_<SAMPLE>.tsv
#   2. Combine per-bin full nulls (bins 6-9) -> data/permutation/counts/trf_bg_signal_<SAMPLE>.tsv
#   3. Re-tag foreground signal files with the NEW bin_ids (interval_id join) —
#      signal quantification itself is bin-independent, so no re-extraction needed.
#   4. Run 3_analyze_period_enrichment.R + period_violin.R
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate r-visualizations
set -euo pipefail

source /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/config_period.sh
init_period_dirs

PER_BIN_DIR="${PERIOD_DATA_DIR}/permutation/per_bin"
SAMPLES=("XG_150" "XG_151" "XG_152" "XG_153")

# ── Step 1: Combine ctrl nulls (bins 1-5, subsample) ──────────────────────────
log "=== Combine control-permutation nulls (bins 1-5) ==="
for s in "${SAMPLES[@]}"; do
  OUT="${COUNTS_DIR}/trf_ctrl_bg_signal_${s}.tsv"
  first=1
  : > "${OUT}.tmp"
  for b in 1 2 3 4 5; do
    f="${PER_BIN_DIR}/ctrl/bin${b}_${s}.tsv"
    [[ -f "$f" ]] || { log "  MISSING ${f}"; continue; }
    if [[ $first -eq 1 ]]; then
      cat "$f" > "${OUT}.tmp"
      first=0
    else
      tail -n +2 "$f" >> "${OUT}.tmp"
    fi
  done
  mv "${OUT}.tmp" "$OUT"
  N=$(tail -n +2 "$OUT" | wc -l)
  log "  [${s}] ${OUT}: ${N} rows"
done

# ── Step 2: Combine full nulls (bins 6-9) ─────────────────────────────────────
log "=== Combine full-set nulls (bins 6-9) ==="
for s in "${SAMPLES[@]}"; do
  OUT="${PERIOD_DATA_DIR}/permutation/counts/trf_bg_signal_${s}.tsv"
  : > "$OUT"
  for b in 6 7 8 9; do
    f="${PER_BIN_DIR}/full/bin${b}_${s}.tsv"
    [[ -f "$f" ]] || { log "  MISSING ${f}"; continue; }
    cat "$f" >> "$OUT"
  done
  log "  [${s}] ${OUT}: $(wc -l < ${OUT}) rows"
done

# ── Step 3: Re-tag foreground signal with new bin_ids ─────────────────────────
log "=== Re-tag foreground signal with new bin_ids ==="
Rscript - "${TRF_BED_CHR}" "${COUNTS_DIR}" <<'REOF'
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(data.table))
bed_all <- args[1]
counts  <- args[2]
samples <- c("XG_150", "XG_151", "XG_152", "XG_153")

bed <- fread(bed_all, header = FALSE,
  col.names = c("chrom","start","end","interval_id","bin_id","period","copies","match"))
map <- bed[, .(interval_id, bin_id)]
message(sprintf("New bin map: %d intervals", nrow(map)))
message(sprintf("Per-bin counts: %s", paste(bed[, .N, by = bin_id][order(bin_id)]$N, collapse = ", ")))

for (s in samples) {
  f <- file.path(counts, paste0("trf_signal_", s, ".tsv"))
  if (!file.exists(f)) { message("  MISSING ", f); next }
  x <- fread(f, header = FALSE,
    col.names = c("chrom","start","end","interval_id","bin_id","period_size",
                  "copies_aligned","match_percent","mean_coverage"))
  stopifnot(!anyDuplicated(x$interval_id), !anyDuplicated(map$interval_id))
  miss <- setdiff(x$interval_id, map$interval_id)
  if (length(miss) > 0) stop(sprintf("  %s: %d interval_ids missing from new BED", s, length(miss)))
  x[, bin_id := NULL]
  x <- merge(x, map, by = "interval_id", sort = FALSE)
  setcolorder(x, c("chrom","start","end","interval_id","bin_id","period_size",
                   "copies_aligned","match_percent","mean_coverage"))
  fwrite(x, f, sep = "\t", col.names = FALSE)
  message(sprintf("  %s: re-tagged %d intervals", s, nrow(x)))
}
REOF

# ── Step 4: Analysis + violin ─────────────────────────────────────────────────
log "=== Step 4: Run analysis R + violin ==="
cd "${PERIOD_DIR}"
Rscript scripts/3_analyze_period_enrichment.R
Rscript scripts/period_violin.R

log "=== Combine + analyze complete ==="
