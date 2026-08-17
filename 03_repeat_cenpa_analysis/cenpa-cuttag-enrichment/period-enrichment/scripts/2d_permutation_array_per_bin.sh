#!/bin/bash
#SBATCH --job-name=trf_perm_array
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --array=1-9
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 2d_permutation_array_per_bin.sh  (sbatch --array=1-9; one job per period bin)
#
# Chromosome- and length-matched permutation null for ONE TRF period bin, on
# MERGED REPEAT ARRAYS (the new unit of observation). Runs in parallel across
# all 9 bins (each bin = one array task).
#
#   Bins 1-5 (>= ~7k merged arrays):  SUBSAMPLE 5,000 arrays × 1,000 iters
#   Bins 6-9 (< 1k merged arrays):    FULL SET        × 1,000 iters
#
# The observed statistic AND each null statistic are computed on MERGED ARRAYS
# (same unit). Shuffles preserve chromosome and merged-array length; shuffled
# arrays may not overlap any TRF locus or assembly gap (same exclusion mask as
# the current pipeline). array_id is carried through the shuffle so the matched
# paired test can join observed signal on the SAME arrays per iteration.
#
# Reads the shared run-length bedGraphs produced by 2e_prepare_bedgraphs.sh
# (already built — reused, not rebuilt).
#
# Output (per bin, per sample):
#   bins 1-5: data/merged/permutation/per_bin/ctrl/bin<N>_<SAMPLE>.tsv  (7 cols WITH header:
#             bin_id, array_id, iter, chrom, start, end, mean_signal)
#   bins 6-9: data/merged/permutation/per_bin/full/bin<N>_<SAMPLE>.tsv  (7 cols no header:
#             chrom, start, end, array_id, bin_id, iter, mean_coverage)
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing
set -euo pipefail

source /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/config_period.sh
init_period_dirs

BIN_ID=$SLURM_ARRAY_TASK_ID
N_ITER=1000
SUBSAMPLE_N=5000
SEED_BASE=20260804

ARRAY_DIR="${PERIOD_DATA_DIR}/merged/arrays"
PERM_DIR="${PERIOD_DATA_DIR}/merged/permutation"
PER_BIN_DIR="${PERM_DIR}/per_bin"
ITER_DIR="${PERM_DIR}/tmp_bin${BIN_ID}"
WORK_DIR="${PERM_DIR}/work_bin${BIN_ID}"
mkdir -p "${PER_BIN_DIR}/ctrl" "${PER_BIN_DIR}/full" "$ITER_DIR" "$WORK_DIR"

CHROM_SIZES="${DATA_DIR}/chrom_sizes_chrOnly.txt"
EXCL_MASK="${PERIOD_DATA_DIR}/permutation/trf_exclusion_mask.bed"
BEDTOOLS_BIN="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/bulk-HiC-processing/bin/bedtools"
BEDGRAPH_DIR="${PERIOD_DATA_DIR}/permutation/bedgraph"

SAMPLES=("XG_150" "XG_151" "XG_152" "XG_153")
SUBSAMPLE_BINS="1 2 3 4 5"

IS_SUBSAMPLE=0
for b in ${SUBSAMPLE_BINS}; do [[ "$b" == "$BIN_ID" ]] && IS_SUBSAMPLE=1; done

log "=== Per-bin permutation on MERGED ARRAYS: bin ${BIN_ID} (subsample=${IS_SUBSAMPLE}) ==="

# ── Step 1: foreground merged arrays for this bin ─────────────────────────────
FG_BIN_BED="${ARRAY_DIR}/bin${BIN_ID}_arrays.bed"
N_FG=$(wc -l < "${FG_BIN_BED}")
if [[ "$N_FG" -eq 0 ]]; then
  log "  Bin ${BIN_ID}: no foreground arrays — writing empty output and exiting"
  for s in "${SAMPLES[@]}"; do
    if [[ "$IS_SUBSAMPLE" -eq 1 ]]; then
      echo -e "bin_id\tarray_id\titer\tchrom\tstart\tend\tmean_signal" > "${PER_BIN_DIR}/ctrl/bin${BIN_ID}_${s}.tsv"
    else
      : > "${PER_BIN_DIR}/full/bin${BIN_ID}_${s}.tsv"
    fi
  done
  rm -rf "$ITER_DIR" "$WORK_DIR"
  exit 0
fi
log "  Foreground merged arrays: ${N_FG}"

# For subsample bins keep the full rows (with array_id) so each iteration can
# draw a FRESH random 5,000-array subsample (reservoir sample, seeded).
FG_FULL_BED="${WORK_DIR}/fg_bin${BIN_ID}_full.bed"
if [[ "$IS_SUBSAMPLE" -eq 1 ]]; then
  cp "${FG_BIN_BED}" "${FG_FULL_BED}"
fi

# ── Step 2: Shuffle N_ITER times (parallel within bin) ───────────────────────
# Subsample bins (1-5): EACH iteration draws a fresh 5,000-array subsample and
#   performs ONE chromosome-/length-matched shuffle per array. The null thus
#   marginalizes over BOTH subsampling and shuffle randomness.
# Full bins (6-9): shuffle the FULL array set each iteration.
log "  Shuffling ${N_ITER} iterations (8 parallel; fresh 5,000-subsample per iter for bins 1-5)..."
shuffle_iter() {
  local iter=$1
  local seed=$((SEED_BASE + BIN_ID*10000 + iter))
  if [[ "$IS_SUBSAMPLE" -eq 1 ]]; then
    awk -v k="$SUBSAMPLE_N" -v s="$seed" 'BEGIN { srand(s) }
      { n++; if (n <= k) pool[n] = $0; else { j = int(rand()*n) + 1; if (j <= k) pool[j] = $0 } }
      END { for (i = 1; i <= k; i++) print pool[i] }' "${FG_FULL_BED}" \
      | "${BEDTOOLS_BIN}" shuffle -i stdin -g "${CHROM_SIZES}" -chrom \
          -excl "${EXCL_MASK}" -seed "${seed}" -noOverlapping 2>/dev/null \
      | awk -v bin="${BIN_ID}" -v iter="${iter}" -v OFS='\t' \
          '{print $1,$2,$3,$4,bin,iter}' \
      > "${ITER_DIR}/iter_${iter}.bed"
  else
    "${BEDTOOLS_BIN}" shuffle -i "${FG_BIN_BED}" -g "${CHROM_SIZES}" -chrom \
      -excl "${EXCL_MASK}" -seed "${seed}" -noOverlapping 2>/dev/null \
      | awk -v bin="${BIN_ID}" -v iter="${iter}" -v OFS='\t' \
          '{print $1,$2,$3,$4,bin,iter}' \
      > "${ITER_DIR}/iter_${iter}.bed"
  fi
}
export -f shuffle_iter
export SEED_BASE BIN_ID IS_SUBSAMPLE SUBSAMPLE_N FG_FULL_BED FG_BIN_BED CHROM_SIZES EXCL_MASK ITER_DIR BEDTOOLS_BIN
seq 1 "${N_ITER}" | xargs -P 8 -n1 -I{} bash -c 'shuffle_iter {}'

BIN_SHUFFLED="${WORK_DIR}/bin${BIN_ID}_all_iters.bed"
cat "${ITER_DIR}"/iter_*.bed > "${BIN_SHUFFLED}"
N_SHUFF=$(wc -l < "${BIN_SHUFFLED}")
log "  Shuffled positions: ${N_SHUFF} (of ${N_ITER}×${N_FG} expected; dropped $((${N_ITER}*N_FG - N_SHUFF)))"

# ── Step 3: Sort + map signal at shuffled positions (all 4 samples) ──────────
sort -k1,1 -k2,2n "${BIN_SHUFFLED}" > "${WORK_DIR}/bin${BIN_ID}_sorted.bed"

for s in "${SAMPLES[@]}"; do
  BG="${BEDGRAPH_DIR}/${s}.bedGraph"
  if [[ ! -f "${BG}" ]]; then
    log "  [${s}] WARNING: bedGraph missing (${BG}) — skipping sample"
    continue
  fi
  MAP_TMP="${WORK_DIR}/map_${s}.tmp"
  "${BEDTOOLS_BIN}" map -a "${WORK_DIR}/bin${BIN_ID}_sorted.bed" -b "${BG}" \
      -c 4 -o mean -null 0 > "${MAP_TMP}"

  if [[ "$IS_SUBSAMPLE" -eq 1 ]]; then
    # shuffled cols: chrom,start,end,array_id,bin_id,iter + mean(col7)
    OUT="${PER_BIN_DIR}/ctrl/bin${BIN_ID}_${s}.tsv"
    echo -e "bin_id\tarray_id\titer\tchrom\tstart\tend\tmean_signal" > "${OUT}"
    awk -v OFS='\t' '{print $5, $4, $6, $1, $2, $3, $7}' "${MAP_TMP}" >> "${OUT}"
  else
    # shuffled cols: chrom,start,end,array_id,bin_id,iter + mean(col7)
    OUT="${PER_BIN_DIR}/full/bin${BIN_ID}_${s}.tsv"
    mv "${MAP_TMP}" "${OUT}"
  fi
  log "  [${s}] wrote $(wc -l < ${OUT}) rows → ${OUT}"
done

# ── Cleanup ──────────────────────────────────────────────────────────────────
rm -rf "$ITER_DIR" "$WORK_DIR"
log "=== Bin ${BIN_ID} complete ==="
