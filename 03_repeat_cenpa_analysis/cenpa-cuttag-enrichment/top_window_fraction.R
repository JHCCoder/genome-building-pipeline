#!/usr/bin/env Rscript
# ============================================================================
# top_window_fraction.R — Per-chromosome CENP-A signal concentration
#
# For each chromosome, compute what fraction of total CENP-A signal falls
# in the single highest 25-kb window. A high fraction indicates one dominant
# centromeric peak; a low fraction suggests distributed or absent signal.
#
# Also reports: top-3-window fraction, top-5-window fraction, and the
# ratio of the top window to the second-highest window (gap statistic).
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

# ----------------------------------------------------------------------------
# Paths and parameters
# ----------------------------------------------------------------------------
BIN_WIDTH <- 25000L
CHROMOSOMES <- paste0("chr", c(1:28, "X", "Y"))

fai_file <- paste0(
  "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/",
  "degu-genome-browser-pythonVersion/",
  "assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta.fai"
)

BW_DIR <- paste0(
  "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/",
  "figure/cenpa-cuttag-centromere/bw_files"
)

bw_files <- list(
  "XG_150" = file.path(BW_DIR, "XG_150.all.bw"),
  "XG_151" = file.path(BW_DIR, "XG_151.all.bw")
)

OUT_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/results"

# ----------------------------------------------------------------------------
# Read chromosome sizes
# ----------------------------------------------------------------------------
chr_sizes <- fread(fai_file, header = FALSE, select = c(1, 2),
                   col.names = c("chrom", "size"))
chr_sizes <- chr_sizes[chrom %in% CHROMOSOMES]
chr_lengths <- setNames(chr_sizes$size, chr_sizes$chrom)

# ----------------------------------------------------------------------------
# Process each BigWig
# ----------------------------------------------------------------------------
all_results <- list()

for (sample_name in names(bw_files)) {

  bw_file <- bw_files[[sample_name]]
  message("Processing ", sample_name, " ...")

  bw <- BigWigFile(bw_file)

  sample_results <- list()

  for (chrom_i in CHROMOSOMES) {

    chr_len <- chr_lengths[[chrom_i]]
    n_bins <- max(1L, as.integer(ceiling(chr_len / BIN_WIDTH)))

    # Summarize BigWig into 25-kb mean bins
    summarized <- summary(
      bw,
      size = n_bins,
      type = "mean",
      which = GRanges(
        seqnames = chrom_i,
        ranges = IRanges(start = 1, end = chr_len)
      )
    )[[1]]

    scores <- as.numeric(summarized$score)
    scores[is.na(scores) | !is.finite(scores)] <- 0
    scores <- pmax(scores, 0)

    if (length(scores) == 0 || sum(scores) == 0) {
      sample_results[[chrom_i]] <- data.table(
        sample = sample_name,
        chrom = chrom_i,
        chr_length_bp = chr_len,
        n_bins = length(scores),
        total_signal = 0,
        top1_score = 0,
        top3_sum = 0,
        top5_sum = 0,
        top1_fraction = NA_real_,
        top3_fraction = NA_real_,
        top5_fraction = NA_real_,
        top1_to_top2_ratio = NA_real_,
        top1_bin_midpoint = NA_real_,
        top1_bin_start = NA_integer_,
        top1_bin_end = NA_integer_
      )
      next
    }

    sorted <- sort(scores, decreasing = TRUE)

    total_signal <- sum(scores)
    top1_score   <- sorted[1]
    top3_sum     <- sum(sorted[1:min(3, length(sorted))])
    top5_sum     <- sum(sorted[1:min(5, length(sorted))])
    top2_score   <- if (length(sorted) >= 2) sorted[2] else 0

    top1_fraction     <- top1_score / total_signal
    top3_fraction     <- top3_sum / total_signal
    top5_fraction     <- top5_sum / total_signal
    top1_to_top2_ratio <- if (top2_score > 0) top1_score / top2_score else Inf

    # Find the genomic position of the top bin
    top1_idx <- which.max(scores)[1]
    bin_start <- summarized@ranges@start[top1_idx]
    bin_end   <- bin_start + summarized@ranges@width[top1_idx] - 1L
    bin_mid   <- floor((bin_start + bin_end) / 2)

    sample_results[[chrom_i]] <- data.table(
      sample = sample_name,
      chrom = chrom_i,
      chr_length_bp = chr_len,
      n_bins = length(scores),
      total_signal = total_signal,
      top1_score = top1_score,
      top3_sum = top3_sum,
      top5_sum = top5_sum,
      top1_fraction = top1_fraction,
      top3_fraction = top3_fraction,
      top5_fraction = top5_fraction,
      top1_to_top2_ratio = top1_to_top2_ratio,
      top1_bin_midpoint = bin_mid,
      top1_bin_start = bin_start,
      top1_bin_end = bin_end
    )
  }

  all_results[[sample_name]] <- rbindlist(sample_results)
}

results <- rbindlist(all_results)

# ----------------------------------------------------------------------------
# Merge across replicates for a consensus metric
# ----------------------------------------------------------------------------
results_wide <- dcast(
  results,
  chrom + chr_length_bp ~ sample,
  value.var = c("top1_fraction", "top3_fraction", "top5_fraction",
                "top1_to_top2_ratio", "top1_bin_midpoint")
)

results_wide[, mean_top1_fraction := (top1_fraction_XG_150 + top1_fraction_XG_151) / 2]
results_wide[, mean_top3_fraction := (top3_fraction_XG_150 + top3_fraction_XG_151) / 2]
results_wide[, mean_top1_to_top2 := (top1_to_top2_ratio_XG_150 + top1_to_top2_ratio_XG_151) / 2]

# Sort by descending top-1 fraction
setorder(results_wide, -mean_top1_fraction)

# ----------------------------------------------------------------------------
# Print summary
# ----------------------------------------------------------------------------
cat("\n")
cat("================================================================================\n")
cat("  Per-chromosome CENP-A signal concentration in the single highest 25-kb window\n")
cat("  Higher fraction = one dominant peak (consistent with single centromere)\n")
cat("================================================================================\n\n")

cat(sprintf("%-8s %12s %12s %12s %12s %10s %12s\n",
            "Chrom", "Top1_frac_1", "Top1_frac_2", "MeanTop1", "Top3_frac",
            "Top1/Top2", "Peak_pos_Mb"))
cat(strrep("-", 82), "\n")

for (r in seq_len(nrow(results_wide))) {
  cat(sprintf("%-8s %12.4f %12.4f %12.4f %12.4f %10.1f %12.2f\n",
              results_wide$chrom[r],
              results_wide$top1_fraction_XG_150[r],
              results_wide$top1_fraction_XG_151[r],
              results_wide$mean_top1_fraction[r],
              results_wide$mean_top3_fraction[r],
              results_wide$mean_top1_to_top2[r],
              results_wide$top1_bin_midpoint_XG_150[r] / 1e6))
}

# ----------------------------------------------------------------------------
# Summary statistics
# ----------------------------------------------------------------------------
cat("\n--- Summary ---\n")
cat(sprintf("Median top-1 fraction (rep1): %.4f\n",
            median(results$top1_fraction[results$sample == "XG_150"], na.rm = TRUE)))
cat(sprintf("Median top-1 fraction (rep2): %.4f\n",
            median(results$top1_fraction[results$sample == "XG_151"], na.rm = TRUE)))
cat(sprintf("Range top-1 fraction: %.4f – %.4f\n",
            min(results_wide$mean_top1_fraction, na.rm = TRUE),
            max(results_wide$mean_top1_fraction, na.rm = TRUE)))
cat(sprintf("Median top-1/top-2 ratio: %.1f\n",
            median(results_wide$mean_top1_to_top2, na.rm = TRUE)))
cat(sprintf("Chromosomes with top-1/top-2 > 2 (dominant peak): %d/%d\n",
            sum(results_wide$mean_top1_to_top2 > 2, na.rm = TRUE),
            nrow(results_wide)))

# Replicate concordance
cor_frac <- cor(results_wide$top1_fraction_XG_150,
                results_wide$top1_fraction_XG_151,
                method = "spearman", use = "complete.obs")
cat(sprintf("Replicate Spearman rho (top-1 fraction): %.3f\n", cor_frac))

# ----------------------------------------------------------------------------
# Write table
# ----------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

fwrite(results,
       file.path(OUT_DIR, "top_window_signal_fraction.csv"))

fwrite(results_wide,
       file.path(OUT_DIR, "top_window_signal_fraction_wide.csv"))

cat("\nResults written to:\n")
cat("  ", file.path(OUT_DIR, "top_window_signal_fraction.csv"), "\n")
cat("  ", file.path(OUT_DIR, "top_window_signal_fraction_wide.csv"), "\n")
cat("\nDone.\n")
