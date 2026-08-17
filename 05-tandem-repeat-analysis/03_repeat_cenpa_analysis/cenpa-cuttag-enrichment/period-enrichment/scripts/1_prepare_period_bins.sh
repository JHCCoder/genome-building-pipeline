#!/bin/bash
#SBATCH --job-name=trf_period_bins
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
# 1_prepare_period_bins.sh
# Extract all TRF repeats from repeat_df_degu.tsv, assign to period bins,
# output sorted BED file with bin labels, and compute bin statistics.
# ============================================================================

source /tscc/nfs/home/jhc103/.bashrc
conda activate r-visualizations

set -euo

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
source "${SCRIPT_DIR}/config_period.sh"
init_period_dirs

log "=== Preparing TRF period-binned BED file ==="

# ============================================================================
# Step 1: Create BED with period bins using R
# ============================================================================
log "Step 1: Reading repeat_df_degu.tsv and assigning period bins"

Rscript --no-save - "${TRF_REPEAT_DF}" "${CHROMOSOMES}" "${TRF_BED_ALL}" "${TRF_BED_CHR}" "${BIN_STATS}" << 'REOF'
args <- commandArgs(trailingOnly = TRUE)
infile       <- args[1]
chromosomes_str <- args[2]
bed_all_out  <- args[3]
bed_chr_out  <- args[4]
stats_out    <- args[5]

chromosomes <- strsplit(chromosomes_str, " ")[[1]]

suppressPackageStartupMessages(library(data.table))

message("Loading TRF repeat dataframe: ", infile)
df <- fread(infile, sep = "\t", header = TRUE, showProgress = TRUE,
            select = c("sequence", "start", "end", "period_size",
                       "copies_aligned", "match_percent", "indel_percent",
                       "alignment_score", "consensus_size"))
message("  Total repeats loaded: ", nrow(df))

# ============================================================================
# Assign period bins
# ============================================================================
# Vectorized bin assignment using findInterval
# Intervals are [break[i], break[i+1]); break values are the LOWER bounds, so
# e.g. breaks c(0,11,51,...) puts period 1-10 -> bin1, 11-50 -> bin2, etc.
# Monomer families 193-195 / 348-349 / 386-390 bp are each isolated.
bin_breaks <- c(0, 11, 51, 193, 196, 348, 350, 386, 391, Inf)
df[, bin_id := findInterval(period_size, bin_breaks, rightmost.closed = TRUE)]

bin_labels <- c(
  "1" = "1-10 bp (microsatellites)",
  "2" = "11-50 bp (minisatellites)",
  "3" = "51-192 bp",
  "4" = "193-195 bp",
  "5" = "196-347 bp",
  "6" = "348-349 bp",
  "7" = "350-385 bp",
  "8" = "386-390 bp",
  "9" = "391+ bp"
)

bin_short_labels <- c(
  "1" = "1-10 bp",
  "2" = "11-50 bp",
  "3" = "51-192 bp",
  "4" = "193-195 bp",
  "5" = "196-347 bp",
  "6" = "348-349 bp",
  "7" = "350-385 bp",
  "8" = "386-390 bp",
  "9" = "391+ bp"
)

n_na <- sum(is.na(df$bin_id))
if (n_na > 0) {
  message("  WARNING: ", n_na, " rows with no bin assigned — removing")
  df <- df[!is.na(bin_id)]
}
message("  After bin assignment: ", nrow(df), " repeats")

# Add bin labels
df[, bin_label := bin_labels[as.character(bin_id)]]
df[, bin_short := bin_short_labels[as.character(bin_id)]]

# ============================================================================
# Compute interval length
# ============================================================================
df[, interval_length := end - start]

# ============================================================================
# Step 2: Convert to BED format (TRF uses 1-based start, BED uses 0-based)
# ============================================================================
df[, bed_start := start - 1L]

# Sort: chr, start
setorder(df, sequence, bed_start)

# Add a unique interval ID
df[, interval_id := .I]

# ============================================================================
# Write full BED (all sequences, including scaffolds)
# BED format: chr, start, end, interval_id, bin_id, period_size, copies, match_pct
# ============================================================================
message("\nWriting full BED (all sequences): ", bed_all_out)
bed_cols <- df[, .(sequence, bed_start, end, interval_id, bin_id,
                    period_size, copies_aligned, match_percent)]
fwrite(bed_cols, bed_all_out, sep = "\t", quote = FALSE, col.names = FALSE)
message("  Done: ", nrow(bed_cols), " intervals")

# ============================================================================
# Write chromosome-only BED
# ============================================================================
message("\nFiltering to main chromosomes for chr-only BED")
df_chr <- df[sequence %in% chromosomes]
message("  Chromosome-only repeats: ", nrow(df_chr))

bed_chr <- df_chr[, .(sequence, bed_start, end, interval_id, bin_id,
                       period_size, copies_aligned, match_percent)]
setorder(bed_chr, sequence, bed_start)
fwrite(bed_chr, bed_chr_out, sep = "\t", quote = FALSE, col.names = FALSE)
message("  Done: ", nrow(bed_chr), " intervals")

# ============================================================================
# Step 3: Compute bin statistics
# ============================================================================
message("\n=== Computing bin statistics ===")

# Overall (all sequences)
stats_all <- df[, .(
  n_loci = .N,
  total_bp = sum(interval_length),
  median_length = as.numeric(median(interval_length)),
  mean_length = as.numeric(mean(interval_length)),
  min_length = min(interval_length),
  max_length = max(interval_length),
  median_copies = as.numeric(median(copies_aligned, na.rm = TRUE)),
  mean_copies = as.numeric(mean(copies_aligned, na.rm = TRUE)),
  median_match = as.numeric(median(match_percent, na.rm = TRUE)),
  mean_match = as.numeric(mean(match_percent, na.rm = TRUE)),
  total_bp_mb = sum(interval_length) / 1e6
), by = .(bin_id)]

stats_all[, bin_label := bin_labels[as.character(bin_id)]]
stats_all[, bin_short := bin_short_labels[as.character(bin_id)]]
stats_all[, category := "all_sequences"]
setorder(stats_all, bin_id)

# Chromosome-only
df_chr_sub <- df[sequence %in% chromosomes]
stats_chr <- df_chr_sub[, .(
  n_loci = .N,
  total_bp = sum(interval_length),
  median_length = as.numeric(median(interval_length)),
  mean_length = as.numeric(mean(interval_length)),
  min_length = min(interval_length),
  max_length = max(interval_length),
  median_copies = as.numeric(median(copies_aligned, na.rm = TRUE)),
  mean_copies = as.numeric(mean(copies_aligned, na.rm = TRUE)),
  median_match = as.numeric(median(match_percent, na.rm = TRUE)),
  mean_match = as.numeric(mean(match_percent, na.rm = TRUE)),
  total_bp_mb = sum(interval_length) / 1e6
), by = .(bin_id)]

stats_chr[, bin_label := bin_labels[as.character(bin_id)]]
stats_chr[, bin_short := bin_short_labels[as.character(bin_id)]]
stats_chr[, category := "chr_only"]
setorder(stats_chr, bin_id)

# Combine
stats_combined <- rbindlist(list(stats_all, stats_chr), use.names = TRUE)
setcolorder(stats_combined, c("category", "bin_id", "bin_short", "bin_label",
                               "n_loci", "total_bp", "total_bp_mb",
                               "median_length", "mean_length", "min_length", "max_length",
                               "median_copies", "mean_copies",
                               "median_match", "mean_match"))

# Write
fwrite(stats_combined, stats_out)
message("Bin statistics saved to: ", stats_out)

# Print summary
cat("\n========================================\n")
cat("Period Bin Statistics (chr only)\n")
cat("========================================\n")
for (i in seq_len(nrow(stats_chr))) {
  s <- stats_chr[i]
  cat(sprintf("\nBin %2s: %s\n", s$bin_id, s$bin_label))
  cat(sprintf("  Loci: %d  |  Total span: %.1f Mb  |  Median length: %.1f kb\n",
      s$n_loci, s$total_bp / 1e6, s$median_length / 1000))
  cat(sprintf("  Mean copies: %.1f  |  Median match: %.1f%%  |  Length range: %d - %d bp\n",
      s$mean_copies, s$median_match, s$min_length, s$max_length))
}
cat("\n========================================\n")

message("\n=== Prepare complete ===")
REOF

log "=== Prepare script complete ==="
log "Full BED: ${TRF_BED_ALL}"
log "Chr-only BED: ${TRF_BED_CHR}"
log "Bin stats: ${BIN_STATS}"
