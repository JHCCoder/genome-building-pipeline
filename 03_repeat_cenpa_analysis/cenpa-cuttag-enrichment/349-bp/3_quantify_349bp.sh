#!/bin/bash
#SBATCH --job-name=349bp_quant
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 3_quantify_349bp.sh
# Count CENP-A/H3K27ac fragments at 349-bp foreground and background intervals
# Uses bedtools coverage (same pattern as V2 pipeline)
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

set -euo

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/349-bp"
source "${SCRIPT_DIR}/config_349bp.sh"
init_349_dirs

log "=== Step 3: Quantify CENP-A signal at 349-bp loci ==="

# Verify inputs
check_file "${FOREGROUND_MERGED}" "foreground merged BED"
check_file "${BG_SHUFFLED}" "shuffled background BED"

# Fragment BEDs
declare -A FRAG_BEDS
FRAG_BEDS["XG_150"]="${FRAGMENT_DIR}/XG_150_fragments.bed"
FRAG_BEDS["XG_151"]="${FRAGMENT_DIR}/XG_151_fragments.bed"
FRAG_BEDS["XG_152"]="${FRAGMENT_DIR}/XG_152_fragments.bed"
FRAG_BEDS["XG_153"]="${FRAGMENT_DIR}/XG_153_fragments.bed"

for sample in "${!FRAG_BEDS[@]}"; do
    check_file "${FRAG_BEDS[$sample]}" "${sample} fragments"
done

# ============================================================================
# Pre-filter fragment BEDs to chr1-28,chrX only (reduces memory for bedtools)
# ============================================================================
log "Filtering fragment BEDs to chr1-28,chrX..."
CHR_PATTERN=$(echo "${CHROMOSOMES_349}" | tr ' ' '|')
FRAG_DIR_349="${DATA_349_DIR}/fragments"
mkdir -p "${FRAG_DIR_349}"

declare -A FRAG_BEDS_FILTERED
for sample in "${!FRAG_BEDS[@]}"; do
    FILTERED="${FRAG_DIR_349}/${sample}_fragments_chrOnly.bed"
    FRAG_BEDS_FILTERED[$sample]="${FILTERED}"
    if [[ -f "${FILTERED}" ]]; then
        log "  ${sample}: already filtered ($(wc -l < ${FILTERED}) lines)"
    else
        grep -E "^(${CHR_PATTERN})" "${FRAG_BEDS[$sample]}" > "${FILTERED}" || true
        log "  ${sample}: filtered $(wc -l < ${FILTERED}) lines (from $(wc -l < ${FRAG_BEDS[$sample]}))"
    fi
done

# ============================================================================
# Count foreground (per chromosome to limit memory)
# ============================================================================
log "Counting foreground intervals (per-chromosome)..."
for sample in "${!FRAG_BEDS_FILTERED[@]}"; do
    OUT_FILE="${COUNTS_349_DIR}/${sample}_349bp_foreground.txt"
    if [[ -f "${OUT_FILE}" ]] && [[ $(wc -l < "${OUT_FILE}") -gt 0 ]]; then
        log "  ${sample} foreground: already exists ($(wc -l < ${OUT_FILE}) lines)"
    else
        log "  ${sample} foreground: counting per chromosome..."
        > "${OUT_FILE}"  # truncate
        for chr in ${CHROMOSOMES_349}; do
            # Extract foreground intervals for this chromosome
            FG_CHR="${DATA_349_DIR}/tmp_fg_${chr}.bed"
            grep "^${chr}[[:space:]]" "${FOREGROUND_MERGED}" > "${FG_CHR}" || true
            if [[ ! -s "${FG_CHR}" ]]; then continue; fi
            # Extract fragments for this chromosome
            FRAG_CHR="${DATA_349_DIR}/tmp_frag_${sample}_${chr}.bed"
            grep "^${chr}[[:space:]]" "${FRAG_BEDS_FILTERED[$sample]}" > "${FRAG_CHR}" || true
            if [[ ! -s "${FRAG_CHR}" ]]; then continue; fi
            # Count for this chromosome
            bedtools coverage -a "${FG_CHR}" -b "${FRAG_CHR}" -counts >> "${OUT_FILE}"
            rm -f "${FG_CHR}" "${FRAG_CHR}"
        done
        log "    Done: $(wc -l < ${OUT_FILE}) intervals"
        rm -f "${DATA_349_DIR}/tmp_fg_"*.bed "${DATA_349_DIR}/tmp_frag_"*.bed
    fi
done

# ============================================================================
# Count background (per chromosome to limit memory)
# ============================================================================
log "Counting background intervals (per-chromosome)..."
for sample in "${!FRAG_BEDS_FILTERED[@]}"; do
    OUT_FILE="${COUNTS_349_DIR}/${sample}_349bp_bg_shuffled.txt"
    if [[ -f "${OUT_FILE}" ]] && [[ $(wc -l < "${OUT_FILE}") -gt 0 ]]; then
        log "  ${sample} background: already exists ($(wc -l < ${OUT_FILE}) lines)"
    else
        log "  ${sample} background: counting per chromosome..."
        > "${OUT_FILE}"
        for chr in ${CHROMOSOMES_349}; do
            BG_CHR="${DATA_349_DIR}/tmp_bg_${chr}.bed"
            grep "^${chr}[[:space:]]" "${BG_SHUFFLED}" > "${BG_CHR}" || true
            if [[ ! -s "${BG_CHR}" ]]; then continue; fi
            FRAG_CHR="${DATA_349_DIR}/tmp_frag_${sample}_${chr}.bed"
            grep "^${chr}[[:space:]]" "${FRAG_BEDS_FILTERED[$sample]}" > "${FRAG_CHR}" || true
            if [[ ! -s "${FRAG_CHR}" ]]; then continue; fi
            bedtools coverage -a "${BG_CHR}" -b "${FRAG_CHR}" -counts >> "${OUT_FILE}"
            rm -f "${BG_CHR}" "${FRAG_CHR}"
        done
        log "    Done: $(wc -l < ${OUT_FILE}) lines"
        rm -f "${DATA_349_DIR}/tmp_bg_"*.bed "${DATA_349_DIR}/tmp_frag_"*.bed
    fi
done

# ============================================================================
# Quick QC
# ============================================================================
log "QC: Fragment count distributions"

conda activate r-visualizations
Rscript --no-save - "${COUNTS_349_DIR}" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
counts_dir <- args[1]

samples <- c("XG_150", "XG_151", "XG_152", "XG_153")
sample_labels <- c("XG_150" = "CENP-A rep1", "XG_151" = "CENP-A rep2",
                   "XG_152" = "H3K27ac", "XG_153" = "H3K27ac rep2")

cat("\n=== Fragment Count Summary ===\n")
for (s in samples) {
    fg_file <- file.path(counts_dir, paste0(s, "_349bp_foreground.txt"))
    bg_file <- file.path(counts_dir, paste0(s, "_349bp_bg_shuffled.txt"))

    if (file.exists(fg_file)) {
        fg <- read.table(fg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        fg_counts <- fg[, ncol(fg)]
        cat(sprintf("\n%s (%s):\n", s, sample_labels[s]))
        cat(sprintf("  Foreground: %d intervals, median fragments = %.1f, mean = %.1f, total = %d\n",
                    length(fg_counts), median(fg_counts), mean(fg_counts), sum(fg_counts)))
        cat(sprintf("              %d / %d intervals with >= 1 fragment\n",
                    sum(fg_counts >= 1), length(fg_counts)))
    }

    if (file.exists(bg_file)) {
        bg <- read.table(bg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        bg_counts <- bg[, ncol(bg)]
        cat(sprintf("  Background: %d intervals, median fragments = %.1f, mean = %.1f\n",
                    length(bg_counts), median(bg_counts), mean(bg_counts)))
    }
}
REOF

log "=== Step 3 complete ==="
log "Count files in: ${COUNTS_349_DIR}/"
ls -lh "${COUNTS_349_DIR}/"
