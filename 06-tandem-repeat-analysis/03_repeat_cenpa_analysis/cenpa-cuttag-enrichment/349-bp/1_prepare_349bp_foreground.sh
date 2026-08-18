#!/bin/bash
#SBATCH --job-name=349bp_prepare
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=you@example.com

# ============================================================================
# 1_prepare_349bp_foreground.sh
# Extract TRF 345-355 bp tandem repeats + centroAnno 349-bp HORs
# Merge into unified 349-bp repeat loci BED file
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

set -euo

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/349-bp"
source "${SCRIPT_DIR}/config_349bp.sh"
init_349_dirs

log "=== Step 1: Prepare 349-bp foreground BED ==="

# ============================================================================
# Part A: Extract TRF 345-355 bp repeats from repeat_df_degu.tsv
# ============================================================================
log "Part A: Extracting TRF repeats with period ${PERIOD_MIN}-${PERIOD_MAX} bp"

if [[ -f "${TRF_349BP_RAW}" ]]; then
    log "TRF 349bp raw BED already exists: ${TRF_349BP_RAW}"
else
    conda activate r-visualizations
    Rscript --no-save - "${TRF_REPEAT_DF}" "${PERIOD_MIN}" "${PERIOD_MAX}" "${CHROMOSOMES_349}" "${TRF_349BP_RAW}" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
pmin <- as.integer(args[2])
pmax <- as.integer(args[3])
chromosomes_str <- args[4]
outfile <- args[5]

chromosomes <- strsplit(chromosomes_str, " ")[[1]]

message("Loading TRF repeat dataframe: ", infile)
df <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message("  Total repeats: ", nrow(df))

# Filter by period and chromosome
df <- df[df$period_size >= pmin & df$period_size <= pmax, ]
message("  After period filter (", pmin, "-", pmax, "): ", nrow(df))

df <- df[df$sequence %in% chromosomes, ]
message("  After chromosome filter: ", nrow(df))

# Convert to BED (TRF uses 1-based start, BED uses 0-based)
bed <- data.frame(
    chrom = df$sequence,
    start = df$start - 1,
    end = df$end,
    stringsAsFactors = FALSE
)
bed <- bed[order(bed$chrom, bed$start), ]

message("Writing TRF 349bp raw BED: ", outfile)
write.table(bed, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("  Done: ", nrow(bed), " intervals")
REOF

    log "TRF 349bp raw BED created: $(wc -l < ${TRF_349BP_RAW}) intervals"
fi

# ============================================================================
# Part B: Extract HORs with ~349 bp period from centroAnno
# ============================================================================
log "Part B: Extracting centroAnno HORs with ~349 bp period"

if [[ -f "${HORS_349BP_RAW}" ]]; then
    log "HORs 349bp raw BED already exists: ${HORS_349BP_RAW}"
else
    # HORs.bed format: chr, HOR_id, start, end, n_monomers, monomer_len, span
    # Column 6 (monomer_len) IS the monomer period in bp — use directly
    conda activate r-visualizations
    Rscript --no-save - "${HORS_BED}" "${CHROMOSOMES_349}" "${PERIOD_MIN}" "${PERIOD_MAX}" "${HORS_349BP_RAW}" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
chromosomes_str <- args[2]
pmin <- as.integer(args[3])
pmax <- as.integer(args[4])
outfile <- args[5]

chromosomes <- strsplit(chromosomes_str, " ")[[1]]

message("Loading HORs: ", infile)
hors <- read.table(infile, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                   col.names = c("chrom", "hor_id", "start", "end", "n_monomers", "monomer_len", "span"))
message("  Total HORs: ", nrow(hors))

# Filter by monomer length (period) — column 6 IS the period
hors <- hors[hors$monomer_len >= pmin & hors$monomer_len <= pmax, ]
message("  After period filter (", pmin, "-", pmax, "): ", nrow(hors))

# Filter by chromosome
hors <- hors[hors$chrom %in% chromosomes, ]
message("  After chromosome filter: ", nrow(hors))

# Write BED (chr, start, end)
bed <- hors[, c("chrom", "start", "end")]
bed <- bed[order(bed$chrom, bed$start), ]

message("Writing HORs 349bp raw BED: ", outfile)
write.table(bed, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message("  Done: ", nrow(bed), " intervals")
REOF

    log "HORs 349bp raw BED created: $(wc -l < ${HORS_349BP_RAW}) intervals"
fi

# ============================================================================
# Part C: Merge all 349-bp intervals
# ============================================================================
log "Part C: Merging TRF + HORs 349-bp intervals (merge dist = ${MERGE_DIST} bp)"

conda activate bulk-HiC-processing

# Concatenate
ALL_RAW="${FOREGROUND_349_DIR}/349bp_all_raw.bed"
cat "${TRF_349BP_RAW}" "${HORS_349BP_RAW}" | sort -k1,1 -k2,2n > "${ALL_RAW}"
log "Combined raw intervals: $(wc -l < ${ALL_RAW})"

# Merge within MERGE_DIST
bedtools merge -i "${ALL_RAW}" -d ${MERGE_DIST} > "${FOREGROUND_MERGED}"
N_MERGED=$(wc -l < "${FOREGROUND_MERGED}")
log "Merged intervals: ${N_MERGED}"

# Add interval IDs
TMP_ID="${FOREGROUND_349_DIR}/349bp_merged_tmp.bed"
awk -v OFS='\t' '{print $1, $2, $3, "r349_" NR, $3-$2}' "${FOREGROUND_MERGED}" > "${TMP_ID}"
mv "${TMP_ID}" "${FOREGROUND_MERGED}"

# ============================================================================
# Part D: QC
# ============================================================================
log "Part D: QC summary"

conda activate r-visualizations
Rscript --no-save - "${FOREGROUND_MERGED}" "${CHROMOSOMES_349}" "${RESULTS_349_DIR}" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
merged_file <- args[1]
chromosomes_str <- args[2]
results_dir <- args[3]

chromosomes <- strsplit(chromosomes_str, " ")[[1]]

df <- read.table(merged_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                 col.names = c("chrom", "start", "end", "id", "size"))

cat("\n=== 349-bp Foreground QC ===\n")
cat("Total merged intervals:", nrow(df), "\n")
cat("Total span:", sum(df$size) / 1e6, "Mb\n")
cat("Median size:", median(df$size) / 1000, "kb\n")
cat("Mean size:", mean(df$size) / 1000, "kb\n")
cat("Size range:", min(df$size), "-", max(df$size), "bp\n\n")

cat("Per chromosome:\n")
for (chr in chromosomes) {
    sub <- df[df$chrom == chr, ]
    if (nrow(sub) > 0) {
        cat(sprintf("  %-6s: %3d intervals, %7.1f kb total\n",
                    chr, nrow(sub), sum(sub$size) / 1000))
    }
}

# Save per-chromosome summary
per_chr <- data.frame(
    chrom = chromosomes,
    n_intervals = sapply(chromosomes, function(c) sum(df$chrom == c)),
    total_bp = sapply(chromosomes, function(c) sum(df$size[df$chrom == c])),
    stringsAsFactors = FALSE
)
write.csv(per_chr, file.path(results_dir, "349bp_per_chromosome.csv"), row.names = FALSE)
cat("\nPer-chromosome summary saved to results/349bp_per_chromosome.csv\n")
REOF

log "=== Step 1 complete ==="
log "Foreground BED: ${FOREGROUND_MERGED} (${N_MERGED} intervals)"
