#!/bin/bash
#SBATCH -J 071826_cenpa_domains
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH --mem=16G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## If parallelization needed:
module load shared

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

# ============================================================================
# 2b_define_domains.sh — Define CENP-A enriched satellite arrays
#
# Strategy:
#   1. Count CENP-A (XG_150/151) fragments at each centroAnno interval
#   2. Normalize CUT&Tag signal and compute mean CENP-A signal per interval
#   3. Flag intervals where mean CENP-A > chromosome median
#   4. Merge flagged intervals at 250 kb → CENP-A enriched satellite arrays
#
# This prevents chr1 from collapsing into one 11.5 Mb super-domain
# (the failure mode of pure distance-based merging).
# ============================================================================

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
source "${SCRIPT_DIR}/config.sh"
init_dirs

THREADS=$SLURM_CPUS_PER_TASK
log "Starting CENP-A-informed domain definition"

SET_A="${FOREGROUND_DIR}/setA_centroAnno_strict.bed"
CHROM_SIZES="${DATA_DIR}/chrom_sizes.txt"
GAP_BED="${EXCLUSION_DIR}/assembly_gaps.bed"

# ============================================================================
# Step 1: Count CENP-A fragments at each centroAnno interval
# ============================================================================
log "=== Step 1: Counting CENP-A fragments at centroAnno intervals ==="

CENPA_SAMPLES=("XG_150" "XG_151")

for sample in "${CENPA_SAMPLES[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    count_file="${DOMAIN_DIR}/${sample}_setA_counts.txt"

    if [[ -f "$count_file" ]] && [[ -s "$count_file" ]]; then
        log "SKIP: $sample Set A counts exist ($(wc -l < $count_file) lines)"
        continue
    fi

    check_file "$frag_bed" "Fragment BED $sample"
    log "Counting $sample fragments at $(wc -l < $SET_A) centroAnno intervals..."

    bedtools coverage -a "$SET_A" -b "$frag_bed" -counts > "$count_file"
    log "  Done: $(wc -l < $count_file) intervals"
done

# ============================================================================
# Step 2: Get library sizes for CUT&Tag signal normalization
# ============================================================================
log "=== Step 2: Computing library sizes ==="

# Count total fragments from fragment BEDs
for sample in "${CENPA_SAMPLES[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    n_frags=$(wc -l < "$frag_bed")
    log "  $sample: $n_frags total fragments"
    echo "$sample $n_frags" > "${DOMAIN_DIR}/${sample}_lib_size.txt"
done

# ============================================================================
# Step 3: Flag intervals by CENP-A signal (R script)
# ============================================================================
log "=== Step 3: Flagging CENP-A-positive intervals ==="

FLAGGED_BED="${DOMAIN_DIR}/cenpa_positive_intervals.bed"
DIAGNOSTIC_CSV="${DOMAIN_DIR}/interval_cenpa_signal.csv"

cat > "${DOMAIN_DIR}/flag_intervals.R" << 'REOF'
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

# Args
args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]

# Read counts
counts_150 <- fread(file.path(base_dir, "XG_150_setA_counts.txt"),
                     header = FALSE, col.names = c("chrom", "start", "end", "id", "count"))
counts_151 <- fread(file.path(base_dir, "XG_151_setA_counts.txt"),
                     header = FALSE, col.names = c("chrom", "start", "end", "id", "count"))

# Read library sizes
lib_150 <- as.numeric(strsplit(readLines(file.path(base_dir, "XG_150_lib_size.txt")), " ")[[1]][2])
lib_151 <- as.numeric(strsplit(readLines(file.path(base_dir, "XG_151_lib_size.txt")), " ")[[1]][2])

# Compute interval length and normalized CUT&Tag signal (use ilen to avoid collision with base::length)
counts_150[, ilen := end - start]
counts_150[, signal_150 := (count / (ilen / 1000)) / (lib_150 / 1e6)]
counts_151[, ilen := end - start]
counts_151[, signal_151 := (count / (ilen / 1000)) / (lib_151 / 1e6)]

# Merge
dt <- merge(counts_150[, .(chrom, start, end, id, ilen, count_150 = count, signal_150)],
            counts_151[, .(chrom, start, end, ilen, count_151 = count, signal_151)],
            by = c("chrom", "start", "end"))

# Pseudocount
all_signal <- c(dt$signal_150, dt$signal_151)
min_nonzero <- min(all_signal[all_signal > 0], na.rm = TRUE)
pseudocount <- min_nonzero / 2

dt[, log2_signal_150 := log2(signal_150 + pseudocount)]
dt[, log2_signal_151 := log2(signal_151 + pseudocount)]
dt[, mean_cenpa := (log2_signal_150 + log2_signal_151) / 2]

# Per-chromosome median CENP-A
chr_medians <- dt[, .(chr_median = median(mean_cenpa, na.rm = TRUE)), by = chrom]

# Flag intervals above chromosome median
dt <- merge(dt, chr_medians, by = "chrom")
dt[, cenpa_positive := mean_cenpa > chr_median]

# Write flagged intervals BED (cenpa_positive = TRUE)
flagged <- dt[(cenpa_positive), .(chrom, start, end, id)]
fwrite(flagged, file.path(base_dir, "cenpa_positive_intervals.bed"),
       sep = "\t", col.names = FALSE)

# Write full diagnostic table
fwrite(dt, file.path(base_dir, "interval_cenpa_signal.csv"))

# Summary
cat(sprintf("Total intervals: %d\n", nrow(dt)))
cat(sprintf("CENP-A-positive intervals: %d (%.1f%%)\n",
            nrow(flagged), 100 * nrow(flagged) / nrow(dt)))
cat(sprintf("Pseudocount: %.6f\n", pseudocount))
cat(sprintf("Mean CENP-A range: %.3f to %.3f\n", min(dt$mean_cenpa), max(dt$mean_cenpa)))

# Per-chromosome summary
chr_summary <- dt[, .(
    n_total = .N,
    n_positive = sum(cenpa_positive),
    pct_positive = round(100 * sum(cenpa_positive) / .N, 1),
    median_cenpa = round(median(mean_cenpa), 3)
), by = chrom][order(chrom)]
cat("\nPer-chromosome flagging:\n")
print(chr_summary, row.names = FALSE)
REOF

conda activate r-visualizations
Rscript "${DOMAIN_DIR}/flag_intervals.R" "${DOMAIN_DIR}"
conda activate bulk-HiC-processing
log "Flagging complete: $(wc -l < $FLAGGED_BED) CENP-A-positive intervals"

# ============================================================================
# Step 4: Merge flagged intervals into domains
# ============================================================================
log "=== Step 4: Merging flagged intervals at ${DOMAIN_MERGE_DIST} bp ==="

MERGED_DOMAINS="${DOMAIN_DIR}/merged_domains_d${DOMAIN_MERGE_DIST}.bed"

# Sort flagged intervals
sort -k1,1V -k2,2n "$FLAGGED_BED" > "${FLAGGED_BED}.sorted"
mv "${FLAGGED_BED}.sorted" "$FLAGGED_BED"

# Merge at 250 kb
bedtools merge -d "$DOMAIN_MERGE_DIST" -i "$FLAGGED_BED" > "${MERGED_DOMAINS}.tmp"

# Add domain IDs and compute metadata
awk -v OFS='\t' '{
    domain_id = sprintf("domain_%04d", NR)
    size = $3 - $2
    print $1, $2, $3, domain_id, size
}' "${MERGED_DOMAINS}.tmp" > "$MERGED_DOMAINS"
rm -f "${MERGED_DOMAINS}.tmp"

n_domains=$(wc -l < "$MERGED_DOMAINS")
log "Merged domains: $n_domains"

# ============================================================================
# Step 5: Count intervals per domain
# ============================================================================
log "=== Step 5: Counting intervals per domain ==="

bedtools intersect -a "$MERGED_DOMAINS" -b "$FLAGGED_BED" -c \
    > "${DOMAIN_DIR}/domains_with_interval_counts.bed"

# Rename columns: chrom, start, end, domain_id, size, n_intervals
awk -v OFS='\t' '{
    print $1, $2, $3, $4, $5, $6
}' "${DOMAIN_DIR}/domains_with_interval_counts.bed" > "${DOMAIN_DIR}/domains_with_interval_counts.bed.tmp"
mv "${DOMAIN_DIR}/domains_with_interval_counts.bed.tmp" "${DOMAIN_DIR}/domains_with_interval_counts.bed"

# Also create a version with all centroAnno intervals (not just flagged)
bedtools intersect -a "$MERGED_DOMAINS" -b "$SET_A" -c \
    > "${DOMAIN_DIR}/domains_with_all_intervals.bed"

# ============================================================================
# Step 6: QC output
# ============================================================================
log "=== Step 6: QC output ==="

# Use R for QC summaries
cat > "${DOMAIN_DIR}/domain_qc.R" << 'REOF'
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]

# Read merged domains
domains <- fread(file.path(base_dir, "domains_with_interval_counts.bed"),
                 header = FALSE, col.names = c("chrom", "start", "end", "domain_id", "size", "n_intervals"))

# Read CENP-A signal (for mean per domain)
signal <- fread(file.path(base_dir, "interval_cenpa_signal.csv"))

# Count flagged intervals per original centroAnno interval
flagged <- fread(file.path(base_dir, "cenpa_positive_intervals.bed"),
                 header = FALSE, col.names = c("chrom", "start", "end", "id"))

cat("=== Domain QC Summary ===\n\n")
cat(sprintf("Total merged domains: %d\n", nrow(domains)))
cat(sprintf("Median domain size: %.1f kb\n", median(domains$size) / 1000))
cat(sprintf("Mean domain size: %.1f kb\n", mean(domains$size) / 1000))
cat(sprintf("Max domain size: %.1f kb\n", max(domains$size) / 1000))
cat(sprintf("Singletons (1 interval): %d (%.1f%%)\n",
            sum(domains$n_intervals == 1),
            100 * sum(domains$n_intervals == 1) / nrow(domains)))
cat(sprintf("Total flagged intervals in domains: %d / %d\n",
            sum(domains$n_intervals), nrow(flagged)))

cat("\nPer-chromosome domain counts:\n")
chr_counts <- domains[, .(n_domains = .N,
                           total_size_mb = round(sum(size) / 1e6, 2),
                           median_size_kb = round(median(size) / 1000, 1),
                           total_intervals = sum(n_intervals)), by = chrom]
chr_counts <- chr_counts[order(chrom)]
print(chr_counts, row.names = FALSE)

# Size distribution
cat("\nDomain size distribution:\n")
size_breaks <- c(0, 10e3, 50e3, 100e3, 250e3, 500e3, 1e6, 5e6, Inf)
size_labels <- c("<10kb", "10-50kb", "50-100kb", "100-250kb", "250-500kb", "500kb-1Mb", "1-5Mb", ">5Mb")
domains[, size_bin := cut(size, breaks = size_breaks, labels = size_labels, right = FALSE)]
size_dist <- domains[, .N, by = size_bin][order(size_bin)]
print(size_dist, row.names = FALSE)

# chr1 check
chr1_domains <- domains[chrom == "chr1"]
cat(sprintf("\nchr1: %d domains (was 1 super-domain with pure distance merge at 500kb)\n", nrow(chr1_domains)))
cat(sprintf("  Largest chr1 domain: %.1f kb, %d intervals\n",
            max(chr1_domains$size) / 1000, chr1_domains[which.max(size), n_intervals]))

# Fraction of flagged intervals per chromosome
cat("\nFraction of centroAnno intervals flagged as CENP-A-positive:\n")
chr_flag <- signal[, .(n_total = .N, n_flagged = sum(cenpa_positive)), by = chrom]
chr_flag[, pct := round(100 * n_flagged / n_total, 1)]
chr_flag <- chr_flag[order(chrom)]
print(chr_flag, row.names = FALSE)

# Write QC table
fwrite(chr_counts, file.path(base_dir, "domains_per_chromosome.csv"))
fwrite(size_dist, file.path(base_dir, "domain_size_distribution.csv"))
REOF

conda activate r-visualizations
Rscript "${DOMAIN_DIR}/domain_qc.R" "${DOMAIN_DIR}" 2>&1 | tee "${DOMAIN_DIR}/domain_qc_summary.txt"
conda activate bulk-HiC-processing

log "Domain QC written to ${DOMAIN_DIR}/domain_qc_summary.txt"

# ============================================================================
# Step 7: Final summary
# ============================================================================
log "=== Domain definition complete ==="
log "Output files:"
log "  Merged domains: $MERGED_DOMAINS ($(wc -l < $MERGED_DOMAINS) domains)"
log "  CENP-A-positive intervals: $FLAGGED_BED ($(wc -l < $FLAGGED_BED) intervals)"
log "  Diagnostic CSV: $DIAGNOSTIC_CSV"
log "  QC summary: ${DOMAIN_DIR}/domain_qc_summary.txt"

log "=== 2b_define_domains.sh DONE ==="
