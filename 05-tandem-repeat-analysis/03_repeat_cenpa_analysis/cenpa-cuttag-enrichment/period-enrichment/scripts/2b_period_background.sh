#!/bin/bash
#SBATCH --job-name=trf_perm_bg
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 2b_period_background.sh
#
# Generates chromosome- and length-matched shuffled background intervals
# for centromeric period bins (6-10), then extracts CENP-A CUT&Tag signal
# at shuffled intervals for empirical P-value computation.
#
# Key design:
#   - Only shuffles bins 6-10 (~3,700 total intervals) — NOT the 1.25M
#     micro/minisatellite intervals (bins 1-5), which serve as negative controls.
#   - 1,000 iterations per bin, chromosome- and length-matched.
#   - Exclusion mask: all TRF intervals + assembly gaps (shuffled intervals
#     cannot land on any known TRF locus or gap).
#   - Signal extracted via bedtools map against per-base coverage bedGraph.
#
# Output: trf_bg_signal_<SAMPLE>.tsv (interval_id, bin_id, iter, mean_coverage)
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

set -euo pipefail

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
source "${SCRIPT_DIR}/config_period.sh"
init_period_dirs

PERM_DIR="${PERIOD_DATA_DIR}/permutation"
BG_DIR="${PERM_DIR}/shuffled_beds"
BG_COUNTS_DIR="${PERM_DIR}/counts"
WORK_DIR="${PERIOD_DATA_DIR}/tmp_bg"
mkdir -p "$BG_DIR" "$BG_COUNTS_DIR" "$WORK_DIR"

FOREGROUND_BED="${TRF_BED_CHR}"
CHROM_SIZES="${DATA_DIR}/chrom_sizes_chrOnly.txt"
GAP_BED="${DATA_DIR}/exclusion/assembly_gaps.bed"

N_ITER=1000
SEED_BASE=20260731
CENTROMERIC_BINS=(6 7 8 9 10)

log "=== Period permutation background ==="
log "Bins to shuffle: ${CENTROMERIC_BINS[*]}"
log "Iterations: ${N_ITER}"
log "Foreground BED: ${FOREGROUND_BED} ($(wc -l < ${FOREGROUND_BED}) intervals)"

# ============================================================================
# Step A: Create exclusion mask (all TRF intervals + assembly gaps)
# ============================================================================
EXCLUSION_MASK="${PERM_DIR}/trf_exclusion_mask.bed"

if [[ ! -f "${EXCLUSION_MASK}" ]]; then
    log "Step A: Creating exclusion mask (TRF intervals + gaps)"

    # TRF intervals: chr, start, end only (first 3 cols)
    cut -f1-3 "${FOREGROUND_BED}" > "${PERM_DIR}/all_trf_intervals.bed"

    # Merge with gaps
    cat "${PERM_DIR}/all_trf_intervals.bed" "${GAP_BED}" \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i stdin \
        > "${EXCLUSION_MASK}"

    N_EXCL=$(wc -l < "${EXCLUSION_MASK}")
    EXCL_BP=$(awk -F'\t' '{sum += $3-$2} END {print sum}' "${EXCLUSION_MASK}")
    log "  Exclusion mask: ${N_EXCL} regions, $(echo "scale=1; $EXCL_BP/1e6" | bc) Mb excluded"
else
    log "Step A: Exclusion mask exists (skip)"
fi

# ============================================================================
# Step B: Generate shuffled intervals per bin
# ============================================================================
BG_CONCAT="${BG_DIR}/all_bins_all_iters.bed"

if [[ ! -f "${BG_CONCAT}" ]]; then
    log "Step B: Generating shuffled backgrounds (${N_ITER} iterations × ${#CENTROMERIC_BINS[@]} bins)"

    for bin_id in "${CENTROMERIC_BINS[@]}"; do
        BIN_BED="${PERM_DIR}/foreground_bin${bin_id}.bed"
        BG_BIN_OUT="${BG_DIR}/bin${bin_id}_all_iters.bed"

        if [[ -f "${BG_BIN_OUT}" ]]; then
            N_BIN=$(wc -l < "${BG_BIN_OUT}")
            log "  Bin ${bin_id}: Already generated (${N_BIN} shuffled intervals)"
            continue
        fi

        # Extract foreground intervals for this bin
        awk -v bin="$bin_id" -F'\t' '$5 == bin' "${FOREGROUND_BED}" \
            | cut -f1-3 > "${BIN_BED}"

        N_FG=$(wc -l < "${BIN_BED}")
        if [[ $N_FG -eq 0 ]]; then
            log "  Bin ${bin_id}: No foreground intervals — skipping"
            continue
        fi
        log "  Bin ${bin_id}: ${N_FG} foreground intervals"

        # Generate shuffled copies
        > "${BG_BIN_OUT}.tmp"  # Initialize empty output
        for iter in $(seq 1 ${N_ITER}); do
            seed=$((SEED_BASE + bin_id * 10000 + iter))

            bedtools shuffle -i "${BIN_BED}" \
                -g "${CHROM_SIZES}" \
                -chrom \
                -excl "${EXCLUSION_MASK}" \
                -seed "${seed}" \
                -noOverlapping 2>/dev/null \
                | awk -v bin="$bin_id" -v iter="$iter" -v OFS='\t' \
                    '{print $1, $2, $3, bin, iter}' \
                >> "${BG_BIN_OUT}.tmp"

            if (( iter % 200 == 0 )); then
                N_DONE=$(wc -l < "${BG_BIN_OUT}.tmp")
                log "    Bin ${bin_id}: iter ${iter}/${N_ITER}, ${N_DONE} total shuffled intervals"
            fi
        done

        mv "${BG_BIN_OUT}.tmp" "${BG_BIN_OUT}"
        N_BIN=$(wc -l < "${BG_BIN_OUT}")
        log "  Bin ${bin_id}: Complete — ${N_BIN} shuffled intervals (${N_ITER} iters × ${N_FG} fg)"
    done

    # Concatenate all bins
    cat "${BG_DIR}"/bin*_all_iters.bed > "${BG_CONCAT}"
    N_TOTAL=$(wc -l < "${BG_CONCAT}")
    log "  Combined background BED: ${N_TOTAL} total shuffled intervals"
    log "  Saved: ${BG_CONCAT}"
else
    N_TOTAL=$(wc -l < "${BG_CONCAT}")
    log "Step B: Background BED exists (${N_TOTAL} intervals) — skip"
fi

# ============================================================================
# Step C: Extract CENP-A signal at shuffled intervals
# ============================================================================
log "Step C: Extracting CENP-A signal at shuffled intervals"

# Sort background BED (required for bedtools map)
BG_SORTED="${WORK_DIR}/bg_sorted.bed"
sort -k1,1 -k2,2n "${BG_CONCAT}" > "${BG_SORTED}"
log "  Background BED sorted"

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    COV_FILE="${COVERAGE_FILES[$i]}"
    OUT_FILE="${BG_COUNTS_DIR}/trf_bg_signal_${SAMPLE}.tsv"

    if [[ -f "${OUT_FILE}" ]]; then
        log "[${SAMPLE}] Background signal already exists: ${OUT_FILE} ($(wc -l < ${OUT_FILE}) rows)"
        continue
    fi

    log "  [${SAMPLE}] Processing..."

    # Convert per-base coverage to bedGraph (same approach as 2_extract_signal.sh)
    BG_COV="${WORK_DIR}/${SAMPLE}_bg_coverage.bedGraph"

    if [[ ! -f "${BG_COV}" ]]; then
        log "    [${SAMPLE}] Converting coverage to bedGraph..."
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
        }' > "${BG_COV}"
        log "    [${SAMPLE}] bedGraph: $(wc -l < ${BG_COV}) blocks"
    fi

    # bedtools map
    log "    [${SAMPLE}] Running bedtools map..."
    bedtools map -a "${BG_SORTED}" -b "${BG_COV}" \
        -c 4 -o mean \
        -null 0 \
        > "${OUT_FILE}.tmp"

    mv "${OUT_FILE}.tmp" "${OUT_FILE}"
    N_BG=$(wc -l < "${OUT_FILE}")
    log "    [${SAMPLE}] Done: ${N_BG} intervals with signal"

    # Clean up bedGraph for this sample
    rm -f "${BG_COV}"
done

# ============================================================================
# Step D: QC summary
# ============================================================================
log "=== Background signal QC ==="

for SAMPLE in "${SAMPLES[@]}"; do
    OUT_FILE="${BG_COUNTS_DIR}/trf_bg_signal_${SAMPLE}.tsv"
    if [[ -f "${OUT_FILE}" ]]; then
        # Columns: chr, start, end, bin_id, iter, mean_coverage
        awk -F'\t' '
        {
            n++; sum+=$6
            if (NR==1) { min=$6; max=$6 }
            if ($6 < min) min=$6
            if ($6 > max) max=$6
            if ($6 > 0) nz++
            bin_count[$4]++
        }
        END {
            printf "  %s: %d intervals, mean_cov=%.4f, nonzero=%d (%.1f%%), range=[%.4f, %.4f]\n",
                   "'${SAMPLE}'", n, sum/(n+0.0001), nz, nz/(n+0.0001)*100, min, max
            printf "    Per bin: "
            for (b in bin_count) printf "bin%s=%d ", b, bin_count[b]
            printf "\n"
        }' "${OUT_FILE}"
    fi
done

# Cleanup
rm -rf "${WORK_DIR}"
log "Cleaned up temporary files"

log "=== Background analysis complete ==="
log "Background counts: ${BG_COUNTS_DIR}/"
log "Shuffled BEDs: ${BG_DIR}/"
