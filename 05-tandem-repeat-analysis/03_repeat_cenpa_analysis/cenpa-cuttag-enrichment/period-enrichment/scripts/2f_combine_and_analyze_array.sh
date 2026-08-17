#!/bin/bash
#SBATCH --job-name=trf_comb_analyze_arr
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
# 2f_combine_and_analyze_array.sh
# Runs AFTER all per-bin MERGED-ARRAY permutation jobs (2d_permutation_array_per_bin.sh,
# array 1-9) and the merged-array signal extraction (2_array_signal.sh) complete.
#
#   1. Combine per-bin ctrl nulls (bins 1-5) -> data/merged/counts/trf_ctrl_bg_signal_array_<SAMPLE>.tsv
#   2. Combine per-bin full nulls (bins 6-9) -> data/merged/permutation/counts/trf_bg_signal_array_<SAMPLE>.tsv
#   3. Run 3_analyze_period_enrichment_array.R (analysis + main figures)
#   4. Run period_violin_array_prepare.R (caches Δ-signal + star positions),
#      then period_violin_array.R (plot-only) for CENP-A and H3K27ac
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate r-visualizations
set -euo pipefail

source /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/config_period.sh
init_period_dirs

PER_BIN_DIR="${PERIOD_DATA_DIR}/merged/permutation/per_bin"
MERGED_COUNTS_DIR="${PERIOD_DATA_DIR}/merged/counts"
FULL_OUT_DIR="${PERIOD_DATA_DIR}/merged/permutation/counts"
mkdir -p "$FULL_OUT_DIR"

SAMPLES=("XG_150" "XG_151" "XG_152" "XG_153")

# ── Step 1: Combine control-permutation nulls (bins 1-5) ─────────────────────
log "=== Combine merged-array ctrl nulls (bins 1-5) ==="
for s in "${SAMPLES[@]}"; do
  OUT="${MERGED_COUNTS_DIR}/trf_ctrl_bg_signal_array_${s}.tsv"
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

# ── Step 2: Combine full-set nulls (bins 6-9) ────────────────────────────────
log "=== Combine merged-array full nulls (bins 6-9) ==="
for s in "${SAMPLES[@]}"; do
  OUT="${FULL_OUT_DIR}/trf_bg_signal_array_${s}.tsv"
  : > "$OUT"
  for b in 6 7 8 9; do
    f="${PER_BIN_DIR}/full/bin${b}_${s}.tsv"
    [[ -f "$f" ]] || { log "  MISSING ${f}"; continue; }
    cat "$f" >> "$OUT"
  done
  log "  [${s}] ${OUT}: $(wc -l < ${OUT}) rows"
done

# ── Step 3: Analysis R (main figures + results) ──────────────────────────────
log "=== Step 3: Run 3_analyze_period_enrichment_array.R ==="
cd "${PERIOD_DIR}"
Rscript scripts/3_analyze_period_enrichment_array.R

# ── Step 4: Violins (CENP-A + H3K27ac) ───────────────────────────────────────
# Prepare caches the per-array Δ-signal and star positions (reads the ~1 GB
# permutation background files); the plot script only reads those caches.
log "=== Step 4: period_violin_array_prepare.R + period_violin_array.R (CENP-A and H3K27ac) ==="
Rscript scripts/period_violin_array_prepare.R CENP-A
Rscript scripts/period_violin_array.R CENP-A
Rscript scripts/period_violin_array_prepare.R H3K27ac
Rscript scripts/period_violin_array.R H3K27ac

log "=== Combine + analyze (merged arrays) complete ==="
