#!/bin/bash
#SBATCH --job-name=trf_signal
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
# 2_extract_signal.sh
# Extract mean CENP-A CUT&Tag signal (per-base coverage) for each TRF interval.
#
# Strategy: Convert per-base coverage to bedGraph (collapse adjacent same-value
# positions) to make bedtools map efficient, then compute mean across each interval.
#
# Input:  trf_chr_only_period_bins.bed (from step 1)
# Output: trf_signal_<SAMPLE>.tsv (interval_id, mean_coverage, n_bases)
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

set -euo pipefail

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
source "${SCRIPT_DIR}/config_period.sh"
init_period_dirs

BED_IN="${TRF_BED_CHR}"
WORK_DIR="${PERIOD_DATA_DIR}/tmp_signal"
mkdir -p "$WORK_DIR"

log "=== Extracting CENP-A signal at TRF intervals ==="
log "Input BED: ${BED_IN} ($(wc -l < ${BED_IN}) intervals)"

# Check that input BED exists and is sorted
if [[ ! -f "${BED_IN}" ]]; then
    log "ERROR: Input BED not found: ${BED_IN}"
    log "Run 1_prepare_period_bins.sh first"
    exit 1
fi

# Check sorting
sort -c -k1,1 -k2,2n "${BED_IN}" 2>/dev/null || {
    log "Sorting BED file..."
    sort -k1,1 -k2,2n "${BED_IN}" > "${BED_IN}.sorted"
    mv "${BED_IN}.sorted" "${BED_IN}"
    log "Sorted."
}

# ============================================================================
# Process each sample
# ============================================================================
for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    COV_FILE="${COVERAGE_FILES[$i]}"
    OUT_FILE="${COUNTS_DIR}/trf_signal_${SAMPLE}.tsv"

    if [[ -f "${OUT_FILE}" ]]; then
        log "[${SAMPLE}] Output already exists: ${OUT_FILE} ($(wc -l < ${OUT_FILE}) rows)"
        continue
    fi

    log "============================================"
    log "[${SAMPLE}] Processing: ${COV_FILE}"
    log "============================================"

    # Step A: Convert per-base coverage to bedGraph (collapse adjacent same values)
    # Format: chr, start (0-based), end, count
    log "[${SAMPLE}] Step A: Converting per-base coverage to bedGraph..."

    BG_FILE="${WORK_DIR}/${SAMPLE}_coverage.bedGraph"

    if [[ ! -f "${BG_FILE}" ]]; then
        # Efficient awk to collapse adjacent same-value positions
        # Coverage format: chr \t pos(1-based) \t count
        # Output bedGraph: chr \t start(0-based) \t end(0-based exclusive) \t count
        # Example: chr1 1 0 -> chr1 0 1 0 (base at 1-based pos 1 = 0-based [0,1))
        zcat "${COV_FILE}" 2>/dev/null | awk '
        BEGIN { FS="\t"; OFS="\t" }
        NR==1 {
            chr=$1; pos=$2; val=$3
            run_start=pos-1; run_end=pos; run_chr=chr; run_val=val
            next
        }
        {
            if ($1 != run_chr || $3 != run_val) {
                # End previous run, start new one
                print run_chr, run_start, run_end, run_val
                run_chr=$1; run_start=$2-1; run_end=$2; run_val=$3
            } else {
                # Extend current run: pos ($2) is 1-based, so exclusive end = pos
                run_end=$2
            }
        }
        END {
            if (NR > 0) print run_chr, run_start, run_end, run_val
        }' > "${BG_FILE}"

        N_BG=$(wc -l < "${BG_FILE}")
        log "[${SAMPLE}] bedGraph blocks: ${N_BG}"
    else
        log "[${SAMPLE}] bedGraph already exists: ${BG_FILE}"
    fi

    # Step B: bedtools map — compute mean coverage per interval
    # -c 4 = column 4 (count value)
    # -o mean = compute mean
    # Also count number of bases covered (for QC)
    log "[${SAMPLE}] Step B: Running bedtools map..."

    bedtools map -a "${BED_IN}" -b "${BG_FILE}" \
        -c 4 -o mean \
        -null 0 \
        > "${OUT_FILE}.tmp"

    # Add header and clean up
    # Output columns from bedtools map:
    # chr, start, end, interval_id, bin_id, period_size, copies, match_pct, mean_coverage

    N_OUT=$(wc -l < "${OUT_FILE}.tmp")
    log "[${SAMPLE}] Mapped intervals: ${N_OUT}"

    # Check for any intervals with no coverage data
    N_NULL=$(awk '$9 == 0 && $4 != ""' "${OUT_FILE}.tmp" | wc -l)
    if [[ ${N_NULL} -gt 0 ]]; then
        log "[${SAMPLE}] WARNING: ${N_NULL} intervals have 0 coverage (no overlapping bedGraph data)"
    fi

    mv "${OUT_FILE}.tmp" "${OUT_FILE}"
    log "[${SAMPLE}] Signal saved to: ${OUT_FILE}"
    log "[${SAMPLE}] Done."
done

# ============================================================================
# QC summary
# ============================================================================
log "=== QC Summary ==="

for SAMPLE in "${SAMPLES[@]}"; do
    OUT_FILE="${COUNTS_DIR}/trf_signal_${SAMPLE}.tsv"
    if [[ -f "${OUT_FILE}" ]]; then
        # Quick stats using awk
        awk -F'\t' '
        NR==1 { next }
        {
            n++; sum+=$9
            if ($9 > 0) nz++
            if (NR==2) { min=$9; max=$9 }
            if ($9 < min) min=$9
            if ($9 > max) max=$9
        }
        END {
            printf "  %s: %d intervals, mean=%.4f, nonzero=%d (%.1f%%), range=[%.4f, %.4f]\n",
                   "'${SAMPLE}'", n, sum/n, nz, nz/n*100, min, max
        }' "${OUT_FILE}"
    fi
done

# Cleanup
rm -rf "${WORK_DIR}"
log "Cleaned up temporary files"

log "=== Signal extraction complete ==="
