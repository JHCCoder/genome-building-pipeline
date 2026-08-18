#!/bin/bash
#SBATCH --job-name=ctrl_perm_bg
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 2c_control_permutation.sh
#
# Permutation null for control bins 1–5.
# Subsamples 1,000 random intervals per bin, then shuffles each interval
# 1,000× (chromosome- and length-matched, excluding TRF regions).
# → 1M null values per bin — directly comparable to bins 6–10 (step 2b).
#
# Output: trf_ctrl_bg_signal_<SAMPLE>.tsv
#   Columns: bin_id, interval_id, iter, chrom, start, end, mean_signal
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate r-visualizations
set -euo pipefail

source /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/config_period.sh
init_period_dirs

N_ITER=1000
SUBSAMPLE_N=1000
SEED_BASE=20260801
CONTROL_BINS=(1 2 3 4 5)

WORK_DIR="${PERIOD_DATA_DIR}/permutation/ctrl_work"
SHUFFLED_DIR="${WORK_DIR}/shuffled"
BG_COUNTS_DIR="${PERIOD_DATA_DIR}/counts"
mkdir -p "$WORK_DIR" "$SHUFFLED_DIR"

CHROM_SIZES="${WORK_DIR}/chrom_sizes.txt"

BEDTOOLS_BIN="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/bulk-HiC-processing/bin/bedtools"

# ── Step 0: Setup ────────────────────────────────────────────────────────────
log "=== Step 0: Setup ==="
log "Bins: ${CONTROL_BINS[*]}"
log "Subsample N: ${SUBSAMPLE_N} intervals per bin"
log "Iterations: ${N_ITER} shuffles per interval"

FAI="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/degu-genome-browser-pythonVersion/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta.fai"
awk '$1 ~ /^chr/ {print $1"\t"$2}' "$FAI" > "$CHROM_SIZES"

# Exclusion mask
EXCL="${WORK_DIR}/exclusion_mask.bed"
if [[ ! -f "$EXCL" ]]; then
  cut -f1-3 "$TRF_BED_CHR" > "$EXCL"
  log "Exclusion mask: $(wc -l < $EXCL) TRF intervals"
fi

# ── Step 1: Subsample 1,000 intervals per bin ─────────────────────────────────
log "=== Step 1: Subsample ${SUBSAMPLE_N} intervals per bin (bins ${CONTROL_BINS[*]}) ==="

SUBSAMPLE_BED="${WORK_DIR}/fg_subsample_ctrl.bed"

if [[ ! -f "${SUBSAMPLE_BED}" ]]; then
  Rscript - "$TRF_BED_CHR" "$SUBSAMPLE_BED" "$SUBSAMPLE_N" "$SEED_BASE" <<'REOF'
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(data.table))
set.seed(as.integer(args[4]))
trf <- fread(args[1], header = FALSE,
  col.names = c("chrom", "start", "end", "interval_id", "bin_id", "period", "copies", "match_pct"))
bins <- 1:5
out <- rbindlist(lapply(bins, function(b) {
  sub <- trf[bin_id == b]
  n_sample <- min(.N, as.integer(args[3]))
  sub[sample(.N, n_sample)]
}))
fwrite(out[, .(chrom, start, end, interval_id, bin_id)], args[2],
       sep = "\t", col.names = FALSE)
for (b in bins) {
  n_b <- sum(out$bin_id == b)
  n_tot <- trf[bin_id == b, .N]
  cat(sprintf("Bin %d: sampled %d / %d intervals\n", b, n_b, n_tot))
}
REOF
  log "Subsampled: $(wc -l < ${SUBSAMPLE_BED}) intervals across ${#CONTROL_BINS[@]} bins"
else
  log "Subsampled BED exists ($(wc -l < ${SUBSAMPLE_BED}) intervals) — skip"
fi

# ── Step 2: Per-bin, 1,000 shuffles per interval ──────────────────────────────
log "=== Step 2: Shuffle ${N_ITER}× per interval ==="

ALL_SHUFFLED="${SHUFFLED_DIR}/all_bins_all_iters.bed"

if [[ ! -f "${ALL_SHUFFLED}" ]]; then
  for bin_id in "${CONTROL_BINS[@]}"; do
    BIN_FG="${WORK_DIR}/fg_bin${bin_id}.bed"
    BIN_OUT="${SHUFFLED_DIR}/bin${bin_id}_all_iters.bed"

    if [[ -f "${BIN_OUT}" ]]; then
      log "  Bin ${bin_id}: exists ($(wc -l < ${BIN_OUT}) rows) — skip"
      continue
    fi

    # Extract subsampled intervals for this bin
    awk -v b="$bin_id" '$5 == b' "$SUBSAMPLE_BED" | cut -f1-4 > "$BIN_FG"
    N_FG=$(wc -l < "$BIN_FG")
    log "  Bin ${bin_id}: ${N_FG} foreground intervals, ${N_ITER} iterations → ~$((N_FG * N_ITER)) null values"

    # Shuffle all intervals in this bin N_ITER times
    > "${BIN_OUT}.tmp"
    for iter in $(seq 1 ${N_ITER}); do
      seed=$((SEED_BASE + bin_id * 10000 + iter))

      ${BEDTOOLS_BIN} shuffle -i "$BIN_FG" \
        -g "$CHROM_SIZES" \
        -chrom \
        -excl "$EXCL" \
        -seed "$seed" \
        -noOverlapping 2>/dev/null \
        | awk -v bin="$bin_id" -v iter="$iter" -v OFS='\t' \
          '{print $1, $2, $3, $4, bin, iter}' \
        >> "${BIN_OUT}.tmp"

      if (( iter % 200 == 0 )); then
        N_DONE=$(wc -l < "${BIN_OUT}.tmp")
        log "    Bin ${bin_id}: iter ${iter}/${N_ITER}, ${N_DONE} total shuffled"
      fi
    done

    mv "${BIN_OUT}.tmp" "${BIN_OUT}"
    log "  Bin ${bin_id}: done — $(wc -l < ${BIN_OUT}) shuffled positions"
  done

  # Concatenate
  cat "${SHUFFLED_DIR}"/bin*_all_iters.bed > "${ALL_SHUFFLED}"
  log "  Combined: $(wc -l < ${ALL_SHUFFLED}) total shuffled intervals"
else
  log "  Combined BED exists ($(wc -l < ${ALL_SHUFFLED}) rows) — skip"
fi

# ── Step 3: Extract signal at shuffled positions ──────────────────────────────
log "=== Step 3: Extract signal at shuffled positions ==="

ALL_SORTED="${WORK_DIR}/all_shuffled_sorted.bed"
sort -k1,1 -k2,2n "${ALL_SHUFFLED}" > "${ALL_SORTED}"
log "  Sorted: $(wc -l < ${ALL_SORTED}) rows"

for i in "${!SAMPLES[@]}"; do
  SAMPLE="${SAMPLES[$i]}"
  COV_FILE="${COVERAGE_FILES[$i]}"
  OUT_FILE="${BG_COUNTS_DIR}/trf_ctrl_bg_signal_${SAMPLE}.tsv"

  if [[ -f "${OUT_FILE}" ]] && [[ $(wc -l < "${OUT_FILE}") -gt 1 ]]; then
    log "  [${SAMPLE}] exists ($(tail -n+2 ${OUT_FILE} | wc -l) rows) — skip"
    continue
  fi

  rm -f "${OUT_FILE}" "${OUT_FILE}.tmp"

  log "  [${SAMPLE}] Converting perbase coverage to bedGraph..."
  BG_COV="${WORK_DIR}/${SAMPLE}_ctrl_coverage.bedGraph"

  # Check for and remove corrupt bedGraph from previous failed run
  if [[ -f "${BG_COV}" ]]; then
    N_BAD=$(awk -F'\t' 'NF!=4 {n++} END {print n+0}' "${BG_COV}")
    if [[ "$N_BAD" -gt 0 ]]; then
      log "    Existing bedGraph has ${N_BAD} corrupt line(s) — removing and regenerating"
      rm -f "${BG_COV}"
    fi
  fi

  if [[ ! -f "${BG_COV}" ]]; then
    zcat "${COV_FILE}" | awk '
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
    }' > "${BG_COV}"

    # Validate: strip any trailing corrupt lines from truncated decompression
    N_BAD=$(awk -F'\t' 'NF!=4 {n++} END {print n+0}' "${BG_COV}")
    if [[ "$N_BAD" -gt 0 ]]; then
      FIRST_BAD=$(awk -F'\t' 'NF!=4 {print NR; exit}' "${BG_COV}")
      log "    WARNING: ${N_BAD} corrupt line(s) at end of bedGraph (zcat may have hit a gzip error)"
      log "    Truncating at line $((FIRST_BAD - 1))"
      head -n $((FIRST_BAD - 1)) "${BG_COV}" > "${BG_COV}.clean"
      mv "${BG_COV}.clean" "${BG_COV}"
    fi
    log "    bedGraph: $(wc -l < ${BG_COV}) blocks"
  else
    log "    bedGraph exists ($(wc -l < ${BG_COV}) blocks) — skip"
  fi

  log "  [${SAMPLE}] ${BEDTOOLS_BIN} map (this may take a while)..."
  ${BEDTOOLS_BIN} map -a "$ALL_SORTED" -b "$BG_COV" -c 4 -o mean -null 0 > "${OUT_FILE}.tmp"

  # Add header, reorder columns
  # shuffled cols: chrom, start, end, interval_id, bin_id, iter
  # map appends: mean_signal (col 7)
  echo -e "bin_id\tinterval_id\titer\tchrom\tstart\tend\tmean_signal" > "${OUT_FILE}"
  awk -v OFS='\t' '{print $5, $4, $6, $1, $2, $3, $7}' "${OUT_FILE}.tmp" >> "${OUT_FILE}"
  rm "${OUT_FILE}.tmp"

  N_ROWS=$(tail -n+2 "${OUT_FILE}" | wc -l)
  log "    Wrote ${N_ROWS} rows → ${OUT_FILE}"
done

# ── Step 4: QC ────────────────────────────────────────────────────────────────
log "=== Step 4: QC ==="
for SAMPLE in "${SAMPLES[@]}"; do
  OUT_FILE="${BG_COUNTS_DIR}/trf_ctrl_bg_signal_${SAMPLE}.tsv"
  if [[ -f "${OUT_FILE}" ]]; then
    N=$(tail -n+2 "${OUT_FILE}" | wc -l)
    if [[ $N -gt 0 ]]; then
      tail -n+2 "${OUT_FILE}" | awk -F'\t' -v s="${SAMPLE}" '
      {
        n++; sum+=$7
        if (NR==1 || $7<min) min=$7
        if ($7>max) max=$7
        if ($7>0) nz++
        bin_ct[$1]++
      }
      END {
        printf "  %s: %d rows, mean=%.4f, nonzero=%d (%.1f%%), range=[%.4f, %.4f]\n",
               s, n, sum/(n+0.0001), nz, nz/(n+0.0001)*100, min, max
        printf "    Per bin: "
        for (b=1; b<=5; b++) printf "bin%d=%d ", b, bin_ct[b]
        printf "\n"
      }'
    else
      log "  [${SAMPLE}]: EMPTY"
    fi
  fi
done

# ── Step 5: Cleanup ───────────────────────────────────────────────────────────
log "=== Step 5: Cleanup ==="
rm -rf "$WORK_DIR"

log "Done. Files: ${BG_COUNTS_DIR}/trf_ctrl_bg_signal_*.tsv"
