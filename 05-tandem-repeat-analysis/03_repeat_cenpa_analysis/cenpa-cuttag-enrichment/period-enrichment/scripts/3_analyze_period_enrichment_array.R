#!/usr/bin/env Rscript
# ============================================================================
# 3_analyze_period_enrichment_array.R
# Period-size stratified CENP-A CUT&Tag enrichment at TRF tandem repeat loci.
# MERGED-ARRAY version: the unit of observation is a merged repeat array
# (bedtools merge -d 0 within each period bin), not an individual TRF interval.
#
# Primary inference: permutation test (observed vs chromosome- and length-
# matched shuffled MERGED ARRAYS, 1,000 iterations). Observed and null
# statistics are both computed on merged arrays (same unit).
# Secondary: Kruskal-Wallis across all bins.
#
# Bin scheme (2026-08-03, 9 bins; monomer families isolated):
#   1 = 1-10 bp     6 = 348-349 bp
#   2 = 11-50 bp    7 = 350-385 bp
#   3 = 51-192 bp   8 = 386-390 bp
#   4 = 193-195 bp  9 = 391+ bp
#   5 = 196-347 bp
# Bins 1-5 (>=7k merged arrays) use a 5,000-array subsample null; bins 6-9 full-set.
#
# All outputs use the "_merged" suffix — the interval-level analysis (no suffix)
# is left untouched for paper writing.
#
# Usage:
#   conda activate r-visualizations
#   Rscript 3_analyze_period_enrichment_array.R
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(patchwork)
  library(ggpubr)
  library(rstatix)
  library(scales)
})

# ============================================================================
# Configuration
# ============================================================================
PERIOD_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PARENT_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
DATA_DIR    <- file.path(PERIOD_DIR, "data")
MERGED_DIR  <- file.path(DATA_DIR, "merged")
ARRAY_DIR   <- file.path(MERGED_DIR, "arrays")
COUNTS_DIR  <- file.path(MERGED_DIR, "counts")
BG_COUNTS_DIR <- file.path(MERGED_DIR, "permutation", "counts")
PLOTS_DIR   <- file.path(PERIOD_DIR, "plots")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_SUFFIX <- "merged"

SAMPLES <- c("XG_150", "XG_151", "XG_152", "XG_153")
SAMPLE_LABELS <- c(
  "XG_150" = "CENP-A rep1",
  "XG_151" = "CENP-A rep2",
  "XG_152" = "H3K27ac",
  "XG_153" = "H3K27ac rep2"
)
SAMPLE_COLORS <- c(
  "CENP-A rep1" = "#2166AC",
  "CENP-A rep2" = "#92C5DE",
  "H3K27ac" = "#B2182B",
  "H3K27ac rep2" = "#D6604D"
)

# Period bin definitions (9 bins; monomer families isolated)
BIN_LABELS_FULL <- c(
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

BIN_SHORT <- c(
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

TESTED_BINS <- 1:9
CENTROMERIC_BINS <- c(6, 7, 8, 9)
CONTROL_BINS <- c(1, 2, 3, 4, 5)

# ── Font size (pt) — single knob; change this one value to scale all text ────
FONT_SIZE <- 10

theme_period <- theme_bw(base_size = FONT_SIZE) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(size = FONT_SIZE - 1, face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    plot.title = element_text(size = FONT_SIZE + 2, face = "bold"),
    plot.subtitle = element_text(size = FONT_SIZE - 1, color = "grey40")
  )

# ============================================================================
# SECTION 0: Redundancy metrics (merged arrays, from 1b)
# ============================================================================
stats_file <- file.path(MERGED_DIR, "merged_array_statistics.csv")
if (file.exists(stats_file)) {
  array_stats <- fread(stats_file)
  message("=== Redundancy metrics (merged arrays, -d 0) ===")
  for (b in TESTED_BINS) {
    r <- array_stats[bin_id == b]
    message(sprintf("  bin %d (%s): %s intervals -> %s arrays (avg %.2f / median %d intervals per array; %s%% in multi arrays)",
                    b, r$bin_label, format(r$n_intervals_total, big.mark=","),
                    format(r$n_arrays, big.mark=","),
                    r$avg_intervals_per_array, r$median_intervals_per_array,
                    sprintf("%.1f", 100 * r$frac_intervals_in_multi)))
  }
} else {
  array_stats <- data.table(bin_id = TESTED_BINS, n_arrays = NA_integer_)
  message("WARNING: merged_array_statistics.csv not found")
}

# ============================================================================
# SECTION 1: Load foreground (merged-array) data
# ============================================================================
message("\n=== SECTION 1: Loading merged-array foreground data ===")

lib_file <- file.path(PARENT_DIR, "data", "qc", "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
} else {
  lib_sizes <- c("XG_150" = 27156788, "XG_151" = 25305681,
                 "XG_152" = 30805412, "XG_153" = 111250740)
}
message("Library sizes (fragments, reused — not recalculated):")
for (s in SAMPLES) message(sprintf("  %s: %d", s, lib_sizes[s]))

signal_list <- list()
for (s in SAMPLES) {
  f <- file.path(COUNTS_DIR, paste0("trf_signal_array_", s, ".tsv"))
  if (!file.exists(f)) {
    message(sprintf("  WARNING: Signal file not found: %s", f))
    next
  }
  x <- fread(f, header = FALSE,
             col.names = c("chrom", "start", "end", "array_id",
                           "bin_id", "n_intervals", "mean_coverage"))
  x[, sample := s]
  x[, array_length := end - start]
  signal_list[[s]] <- x
  message(sprintf("  %s: %d merged arrays", s, nrow(x)))
}
signal <- rbindlist(signal_list)
message(sprintf("  Total foreground: %d rows (merged arrays)", nrow(signal)))

if (nrow(signal) == 0) {
  stop("No merged-array foreground signal data. Run 2_array_signal.sh first.")
}

signal[, sample_label := factor(SAMPLE_LABELS[sample],
                                 levels = c("CENP-A rep1", "CENP-A rep2",
                                            "H3K27ac", "H3K27ac rep2"))]
signal[, bin_factor := factor(bin_id, levels = 1:9)]

# ============================================================================
# SECTION 2: Normalize
# ============================================================================
message("\n=== SECTION 2: Normalizing signal (merged arrays) ===")

signal[, lib_size := lib_sizes[sample]]
signal[, norm_signal := mean_coverage / (lib_size / 1e6)]

pseudo <- min(signal$norm_signal[signal$norm_signal > 0], na.rm = TRUE) / 2
message(sprintf("Pseudocount (merged arrays): %.4e", pseudo))

signal[, log2_signal := log2(norm_signal + pseudo)]

# ============================================================================
# SECTION 3: Load background (permuted) data — merged arrays
# ============================================================================
message("\n=== SECTION 3: Loading background (permuted) data ===")

bg_list <- list()
has_bg <- FALSE
for (s in SAMPLES) {
  f <- file.path(BG_COUNTS_DIR, paste0("trf_bg_signal_array_", s, ".tsv"))
  if (!file.exists(f)) {
    message(sprintf("  WARNING: No full-set background data for %s", s))
    next
  }
  # Columns: chrom, start, end, array_id, bin_id, iter, mean_coverage (7 cols, no header)
  x <- fread(f, header = FALSE,
             col.names = c("chrom", "start", "end", "array_id", "bin_id", "iter", "mean_coverage"))
  x[, sample := s]
  x[, norm_signal := mean_coverage / (lib_sizes[s] / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
  n_bg <- nrow(x)
  n_bins <- uniqueN(x$bin_id)
  bg_list[[s]] <- x[, .(sample, chrom, bin_id, iter, log2_signal)]
  has_bg <- TRUE
  message(sprintf("  %s: %d full-set background rows across %d bins", s, n_bg, n_bins))
}

ctrl_has_data <- FALSE
for (s in SAMPLES) {
  f <- file.path(COUNTS_DIR, paste0("trf_ctrl_bg_signal_array_", s, ".tsv"))
  if (!file.exists(f)) {
    message(sprintf("  NOTE: No control-permutation background for %s", s))
    next
  }
  # Columns (7 cols WITH header): bin_id, array_id, iter, chrom, start, end, mean_signal
  x <- fread(f, header = TRUE)
  setnames(x, "mean_signal", "mean_coverage")
  x[, sample := s]
  x[, array_length := end - start]
  # Join OBSERVED signal for the same array_ids, so the matched permutation test
  # compares observed vs shuffled medians computed on the SAME merged arrays.
  obs <- signal[sample == s, .(array_id, obs_coverage = mean_coverage)]
  x <- merge(x, obs, by = "array_id", all.x = TRUE)
  n_miss <- sum(is.na(x$obs_coverage))
  if (n_miss > 0) message(sprintf("    WARNING: %d/%d ctrl rows missing observed coverage", n_miss, nrow(x)))
  x[, norm_obs := obs_coverage / (lib_sizes[s] / 1e6)]
  x[, log2_obs := log2(norm_obs + pseudo)]
  x[, norm_signal := mean_coverage / (lib_sizes[s] / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
  x_compact <- x[, .(sample, chrom, bin_id, iter, log2_signal, log2_obs)]
  if (is.null(bg_list[[s]])) {
    bg_list[[s]] <- x_compact
  } else {
    bg_list[[s]] <- rbindlist(list(bg_list[[s]], x_compact), fill = TRUE)
  }
  has_bg <- TRUE
  ctrl_has_data <- TRUE
  message(sprintf("  %s: %d control-permutation rows across %d bins", s, nrow(x), uniqueN(x$bin_id)))
}

if (has_bg) {
  bg <- rbindlist(bg_list, fill = TRUE)
  message(sprintf("  Total background: %d rows", nrow(bg)))
  message("  (compact: sample, chrom, bin_id, iter, log2_signal, log2_obs)")
  bg_iters <- bg[, .(n_iters = uniqueN(iter)), by = .(sample, bin_id)]
  message(sprintf("  Iterations per bin per sample: %d - %d", min(bg_iters$n_iters), max(bg_iters$n_iters)))
}

# ============================================================================
# SECTION 4: Permutation test (primary inference) — merged arrays
# ============================================================================
message("\n=== SECTION 4: Permutation test (merged arrays) ===")

if (has_bg) {
  fg_bin <- signal[bin_id %in% TESTED_BINS,
                    .(fg_median = median(log2_signal, na.rm = TRUE),
                      fg_mean = mean(log2_signal, na.rm = TRUE),
                      n_fg = .N),
                    by = .(sample, bin_id)]

  # Matched (paired) permutation: per-iteration observed AND shuffled medians
  # computed on the SAME merged arrays.
  #   bins 1-5: obs_median = median observed log2 signal over the 5,000 arrays
  #             subsampled AND shuffled in that iteration (array_id join).
  #   bins 6-9: shuffled set = full array set; observed = full-bin median (constant).
  bg_bin_iter <- bg[bin_id %in% TESTED_BINS,
                     .(bg_median = median(log2_signal, na.rm = TRUE),
                       bg_mean = mean(log2_signal, na.rm = TRUE),
                       obs_median = median(log2_obs, na.rm = TRUE)),
                     by = .(sample, bin_id, iter)]
  bg_bin_iter <- merge(bg_bin_iter, fg_bin[, .(sample, bin_id, fg_median)],
                       by = c("sample", "bin_id"))
  bg_bin_iter[, obs_med := ifelse(is.na(obs_median), fg_median, obs_median)]
  bg_bin_iter[, D := obs_med - bg_median]

  perm_pvals <- bg_bin_iter[, .(
    fg_median = fg_median[1],
    p_enrich = mean(D <= 0),
    p_deplete = mean(D >= 0),
    n_iters = .N,
    bg_mean_of_medians = mean(bg_median),
    bg_sd_of_medians = sd(bg_median),
    obs_mean_of_medians = mean(obs_med),
    paired_effect = mean(D),
    paired_effect_sd = sd(D),
    z_score = mean(D) / (sd(D) + 1e-10)
  ), by = .(sample, bin_id)]
  perm_pvals[paired_effect_sd == 0, z_score := NA_real_]

  perm_pvals[, bin_label := BIN_LABELS_FULL[as.character(bin_id)]]
  perm_pvals[, sample_label := SAMPLE_LABELS[sample]]
  perm_pvals[, p_two_sided := pmin(p_enrich, p_deplete) * 2]
  perm_pvals[, p_two_sided := pmin(p_two_sided, 1.0)]
  perm_pvals[, p_adj := p.adjust(p_two_sided, method = "BH"), by = sample]
  perm_pvals[, sig := ifelse(p_adj < 0.001, "***",
                      ifelse(p_adj < 0.01, "**",
                      ifelse(p_adj < 0.05, "*",
                      ifelse(p_adj < 0.1, ".", "ns"))))]

  # Add true n_arrays per bin (from the merge statistics)
  if ("n_arrays" %in% names(array_stats)) {
    perm_pvals <- merge(perm_pvals, array_stats[, .(bin_id, n_arrays, n_intervals_total)],
                        by = "bin_id", all.x = TRUE)
  }

  message(sprintf("\nperm_pvals dimensions: %d rows x %d cols", nrow(perm_pvals), ncol(perm_pvals)))

  message("\nPermutation test results (empirical P, one-sided enrichment; merged arrays):")
  message("===============================================================")
  for (s in c("XG_150", "XG_151")) {
    message(sprintf("\n%s:", SAMPLE_LABELS[s]))
    sub <- perm_pvals[sample == s][order(bin_id)]
    if (nrow(sub) == 0) next
    for (i in seq_len(nrow(sub))) {
      r <- sub[i]
      zs <- if (is.na(r$z_score)) "  —" else sprintf("%.2f", r$z_score)
      message(sprintf("  Bin %s (%s): fg_median=%.3f, bg_null=%.3f±%.3f, D=%.3f±%.3f, z=%s, P_enrich=%.3f, P_adj=%.3f %s",
                      r$bin_id, r$bin_label, r$fg_median,
                      r$bg_mean_of_medians, r$bg_sd_of_medians,
                      r$paired_effect, r$paired_effect_sd,
                      zs, r$p_enrich, r$p_adj, r$sig))
    }
  }
} else {
  message("  SKIP: No background data available.")
}

# ============================================================================
# SECTION 5: Kruskal-Wallis (secondary)
# ============================================================================
message("\n=== SECTION 5: Kruskal-Wallis (merged arrays) ===")

kw_results <- list()
for (s in c("XG_150", "XG_151")) {
  sub <- signal[sample == s & bin_id <= 9]
  kw <- kruskal.test(log2_signal ~ bin_factor, data = sub)
  message(sprintf("%s: chi2=%.2f, df=%d, P=%.3e",
                  SAMPLE_LABELS[s], kw$statistic, kw$parameter, kw$p.value))
  dunn_res <- sub %>%
    dunn_test(log2_signal ~ bin_factor, p.adjust.method = "bonferroni") %>%
    as.data.table()
  kw_results[[s]] <- list(kw = kw, dunn = dunn_res)
}

# ============================================================================
# SECTION 6: Main figure — period-stratified enrichment (merged arrays)
# ============================================================================
message("\n=== SECTION 6: Main figure ===")

if (has_bg) {
  pval_cols <- intersect(c("sample", "bin_id", "sig", "p_adj"), names(perm_pvals))
  if (length(pval_cols) >= 4 && nrow(perm_pvals) > 0) {
    fg_bin_plot <- merge(fg_bin, perm_pvals[, .SD, .SDcols = pval_cols],
                         by = c("sample", "bin_id"), all.x = TRUE)
  } else {
    fg_bin_plot <- copy(fg_bin)
    fg_bin_plot[, sig := " "]
    fg_bin_plot[, p_adj := NA_real_]
  }
  if (!"sig" %in% names(fg_bin_plot)) fg_bin_plot[, sig := " "]
  fg_bin_plot[is.na(sig), sig := " "]
  fg_bin_plot[, sample_label := factor(SAMPLE_LABELS[sample],
                                        levels = c("CENP-A rep1", "CENP-A rep2"))]

  bg_bin_dist <- bg_bin_iter[, .(
    median_of_median = median(bg_median),
    q025 = quantile(bg_median, 0.025),
    q975 = quantile(bg_median, 0.975),
    q25 = quantile(bg_median, 0.25),
    q75 = quantile(bg_median, 0.75)
  ), by = .(sample, bin_id)]
  bg_bin_dist[, sample_label := factor(SAMPLE_LABELS[sample],
                                        levels = c("CENP-A rep1", "CENP-A rep2"))]

  p_perm <- ggplot() +
    geom_rect(data = bg_bin_dist,
              aes(xmin = as.numeric(as.character(bin_id)) - 0.35,
                  xmax = as.numeric(as.character(bin_id)) + 0.35,
                  ymin = q025, ymax = q975),
              fill = "grey85", alpha = 0.7) +
    geom_rect(data = bg_bin_dist,
              aes(xmin = as.numeric(as.character(bin_id)) - 0.35,
                  xmax = as.numeric(as.character(bin_id)) + 0.35,
                  ymin = q25, ymax = q75),
              fill = "grey70", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(data = fg_bin_plot,
               aes(x = as.numeric(as.character(bin_id)), y = fg_median,
                   color = sample_label),
               size = 3, shape = 18) +
    geom_text(data = fg_bin_plot[sig != " "],
              aes(x = as.numeric(as.character(bin_id)),
                  y = fg_median + 0.2, label = sig),
              size = 5, color = "red", fontface = "bold") +
    geom_point(data = bg_bin_dist,
               aes(x = as.numeric(as.character(bin_id)), y = median_of_median),
               shape = 3, size = 1.5, color = "grey40") +
    facet_wrap(~ sample_label, ncol = 1) +
    scale_color_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_x_continuous(breaks = TESTED_BINS, labels = BIN_SHORT[as.character(TESTED_BINS)]) +
    labs(
      x = "TRF repeat period bin",
      y = expression("Median log"[2] * "(normalized CUT&Tag signal)")
    ) +
    theme_period

  p_perm_full <- p_perm + theme(axis.text.x = element_text(size = FONT_SIZE - 2, angle = 0, hjust = 0.5))
} else {
  message("No background data — using boxplot-only figure")
  cenpa_signal <- signal[sample %in% c("XG_150", "XG_151") & bin_id <= 9]
  p_perm_full <- ggplot(cenpa_signal, aes(x = bin_factor, y = log2_signal)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_boxplot(aes(fill = bin_id %in% TESTED_BINS),
                 outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.25) +
    facet_wrap(~ sample_label, ncol = 1) +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "grey80"), guide = "none") +
    scale_x_discrete(labels = BIN_SHORT[as.character(1:9)]) +
    labs(x = "TRF repeat period", y = expression(log[2] * "(normalized signal)")) +
    theme_period +
    theme(axis.text.x = element_text(size = FONT_SIZE - 3))
}

ggsave(file.path(PLOTS_DIR, paste0("period_permutation_enrichment_", OUT_SUFFIX, ".pdf")),
       p_perm_full, width = 8, height = 7)
ggsave(file.path(PLOTS_DIR, paste0("period_permutation_enrichment_", OUT_SUFFIX, ".png")),
       p_perm_full, width = 8, height = 7, dpi = 300)
message("Saved: period_permutation_enrichment_merged.pdf/png")

# Panel B: All bins overview (foreground only)
fg_all_bins <- signal[bin_id <= 9, .(
  median_log2 = median(log2_signal, na.rm = TRUE),
  mean_log2 = mean(log2_signal, na.rm = TRUE),
  se_log2 = sd(log2_signal, na.rm = TRUE) / sqrt(.N),
  n = .N
), by = .(sample, bin_id)]
fg_all_bins[, sample_label := SAMPLE_LABELS[sample]]
fg_all_bins[, bin_id_num := as.numeric(as.character(bin_id))]

p_overview <- ggplot(fg_all_bins, aes(x = bin_id_num, y = mean_log2,
                                       color = sample_label, group = sample_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 5.5, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  annotate("text", x = 1, y = Inf, label = "Micro/mini/satellite bins", vjust = 1.5, hjust = 0,
           size = 3, color = "grey50", fontface = "italic") +
  annotate("text", x = 7, y = Inf, label = "Centromeric monomer bins", vjust = 1.5, hjust = 0,
           size = 3, color = "grey50", fontface = "italic") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n), alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_log2 - se_log2, ymax = mean_log2 + se_log2),
                width = 0.3, linewidth = 0.5) +
  scale_color_manual(values = SAMPLE_COLORS, name = NULL) +
  scale_x_continuous(breaks = 1:9, labels = BIN_SHORT[as.character(1:9)]) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(x = "TRF repeat period", y = expression("Mean log"[2] * "(normalized signal) ± SE")) +
  theme_period +
  theme(axis.text.x = element_text(size = 7))

ggsave(file.path(PLOTS_DIR, paste0("period_all_bins_overview_", OUT_SUFFIX, ".pdf")),
       p_overview, width = 10, height = 5)
ggsave(file.path(PLOTS_DIR, paste0("period_all_bins_overview_", OUT_SUFFIX, ".png")),
       p_overview, width = 10, height = 5, dpi = 300)
message("Saved: period_all_bins_overview_merged.pdf/png")

# ============================================================================
# SECTION 7: Supplementary figures
# ============================================================================
message("\n=== SECTION 7: Supplementary figures ===")

# Supp A: Z-score heatmap per chromosome per bin (CENP-A rep1)
# NOTE: the PAPER version of this heatmap is produced by
# scripts/supp_permutation_zscore_heatmap_merged.R (bins 3-9 only,
# RED=enriched / BLUE=depleted, per-chromosome permutation significance stars,
# chr1 at top, x-axis at bottom). It is intentionally NOT written here so this
# analysis cannot overwrite that file.
message("Supp A heatmap: produced by scripts/supp_permutation_zscore_heatmap_merged.R (not written here)")

# Supp B: 348-349 bp monomer (bin 6) null distribution detail (merged arrays)
if (has_bg && nrow(perm_pvals) > 0 && "sig" %in% names(perm_pvals)) {
  bg_bin6 <- bg_bin_iter[bin_id == 6]
  fg_bin6 <- fg_bin[bin_id == 6]
  fg_bin6[, sample_label := factor(SAMPLE_LABELS[sample],
                                    levels = c("CENP-A rep1", "CENP-A rep2",
                                               "H3K27ac", "H3K27ac rep2"))]
  bg_bin6[, sample_label := factor(SAMPLE_LABELS[sample],
                                    levels = c("CENP-A rep1", "CENP-A rep2",
                                               "H3K27ac", "H3K27ac rep2"))]

  p_null <- ggplot(bg_bin6[sample %in% c("XG_150", "XG_151")]) +
    geom_histogram(aes(x = bg_median, fill = sample_label),
                   bins = 60, alpha = 0.6, position = "identity") +
    geom_vline(data = fg_bin6[sample %in% c("XG_150", "XG_151")],
               aes(xintercept = fg_median, color = sample_label),
               linewidth = 1.2) +
    geom_text(data = perm_pvals[bin_id == 6 & sample %in% c("XG_150", "XG_151")],
              aes(x = Inf, y = Inf,
                  label = sprintf("P_enrich = %.3f\nP_adj = %.3f %s",
                                  p_enrich, p_adj, sig)),
              hjust = 1.05, vjust = 1.1, size = 3, fontface = "bold") +
    facet_wrap(~ sample_label, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_color_manual(values = SAMPLE_COLORS, guide = "none") +
    labs(x = expression("Median log"[2] * "(normalized signal) per iteration"),
         y = "Frequency (1,000 iterations)") +
    theme_period

  ggsave(file.path(PLOTS_DIR, paste0("supp_bin6_null_distribution_", OUT_SUFFIX, ".pdf")),
         p_null, width = 7, height = 6)
  ggsave(file.path(PLOTS_DIR, paste0("supp_bin6_null_distribution_", OUT_SUFFIX, ".png")),
         p_null, width = 7, height = 6, dpi = 300)
  message("Saved: supp_bin6_null_distribution_merged.pdf/png")
}

# Supp C: Copy number vs signal (within 348-349 bp monomer, bin 6) — per ARRAY
TRF_BED <- file.path(DATA_DIR, "trf_chr_only_period_bins.bed")
if (file.exists(TRF_BED)) {
  trf_all <- fread(TRF_BED, header = FALSE,
    col.names = c("chrom", "start", "end", "interval_id", "bin_id",
                  "period", "copies", "match_pct"))
  merged_bed <- fread(file.path(ARRAY_DIR, "merged_arrays_all_bins.bed"), header = FALSE,
    col.names = c("chrom", "start", "end", "array_id", "bin_id", "n_intervals"))
  # Per-array copy number / match% = aggregate over the array's constituent TRF
  # intervals (foverlaps; arrays are disjoint, so each interval maps to 1 array).
  trf6 <- trf_all[bin_id == 6, .(chrom, start, end, interval_id, copies_aligned = copies, match_percent = match_pct)]
  merged6 <- merged_bed[bin_id == 6, .(chrom, start, end, array_id)]
  setkey(trf6, chrom, start, end)
  setkey(merged6, chrom, start, end)
  mem <- foverlaps(trf6, merged6, type = "any")
  arr6_copies <- mem[, .(copies_total = sum(copies_aligned, na.rm = TRUE),
                         match_mean = mean(match_percent, na.rm = TRUE),
                         n_int = .N), by = array_id]
  bin6_signal <- signal[bin_id == 6 & sample %in% c("XG_150", "XG_151")]
  bin6_signal <- merge(bin6_signal, arr6_copies, by = "array_id", all.x = TRUE)
  bin6_signal[, sample_label := factor(SAMPLE_LABELS[sample],
                                        levels = c("CENP-A rep1", "CENP-A rep2"))]

  p_copies <- ggplot(bin6_signal, aes(x = log10(copies_total + 1), y = log2_signal)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_point(aes(color = chrom), size = 0.5, alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "#2166AC", linewidth = 0.8) +
    facet_wrap(~ sample_label, nrow = 1) +
    scale_color_discrete(guide = "none") +
    labs(x = expression(log[10] * "(total array copies + 1)"),
         y = expression(log[2] * "(normalized signal)")) +
    theme_period

  ggsave(file.path(PLOTS_DIR, paste0("supp_349bp_copies_vs_signal_", OUT_SUFFIX, ".pdf")),
         p_copies, width = 9, height = 4)
  ggsave(file.path(PLOTS_DIR, paste0("supp_349bp_copies_vs_signal_", OUT_SUFFIX, ".png")),
         p_copies, width = 9, height = 4, dpi = 300)
  message("Saved: supp_349bp_copies_vs_signal_merged.pdf/png")

  # Supp D: Match% (repeat homogeneity) vs signal — per array
  p_match <- ggplot(bin6_signal, aes(x = match_mean, y = log2_signal)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_point(aes(color = chrom), size = 0.5, alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "#2166AC", linewidth = 0.8) +
    facet_wrap(~ sample_label, nrow = 1) +
    scale_color_discrete(guide = "none") +
    labs(x = "Mean match percent (%)", y = expression(log[2] * "(normalized signal)")) +
    theme_period

  ggsave(file.path(PLOTS_DIR, paste0("supp_349bp_match_vs_signal_", OUT_SUFFIX, ".pdf")),
         p_match, width = 9, height = 4)
  ggsave(file.path(PLOTS_DIR, paste0("supp_349bp_match_vs_signal_", OUT_SUFFIX, ".png")),
         p_match, width = 9, height = 4, dpi = 300)
  message("Saved: supp_349bp_match_vs_signal_merged.pdf/png")
}

# ============================================================================
# SECTION 7b: Zero-inflation & non-independence supplement (merged arrays)
# ============================================================================
message("\n=== SECTION 7b: Zero-inflation supplement (merged arrays) ===")

if (has_bg) {
  null_base <- bg_bin_iter[, .(null_median = mean(bg_median)), by = .(sample, bin_id)]
  fg_all <- signal[bin_id %in% TESTED_BINS]
  fg_all <- merge(fg_all, null_base, by = c("sample", "bin_id"))
  fg_all[, delta := log2_signal - null_median]

  null_p95 <- bg[bin_id %in% TESTED_BINS,
                 .(null_p95 = as.numeric(quantile(log2_signal, 0.95, na.rm = TRUE))),
                 by = .(sample, bin_id)]
  fg_all <- merge(fg_all, null_p95, by = c("sample", "bin_id"))

  zmeas <- fg_all[, .(
    n = .N,
    frac_nonzero = mean(norm_signal > 0, na.rm = TRUE),
    frac_delta_gt0 = mean(delta > 0, na.rm = TRUE),
    p75_delta = as.numeric(quantile(delta, 0.75, na.rm = TRUE)),
    p90_delta = as.numeric(quantile(delta, 0.90, na.rm = TRUE)),
    frac_fg_above_null_p95 = mean(log2_signal > null_p95, na.rm = TRUE)
  ), by = .(sample, bin_id)]
  # TRUE merged-array counts (from 1b statistics — replaces the old buggy
  # uniqueN-over-per-chromosome-ID computation)
  if ("n_arrays" %in% names(array_stats)) {
    zmeas <- merge(zmeas, array_stats[, .(bin_id, n_arrays, n_intervals_total)],
                   by = "bin_id", all.x = TRUE)
  } else {
    zmeas[, `:=`(n_arrays = NA_integer_, n_intervals_total = NA_integer_)]
  }
  zmeas <- merge(zmeas, null_p95, by = c("sample", "bin_id"))
  zmeas[, bin_label := BIN_LABELS_FULL[as.character(bin_id)]]
  zmeas[, sample_label := SAMPLE_LABELS[sample]]
  setcolorder(zmeas, c("sample", "sample_label", "bin_id", "bin_label", "n",
                       "n_intervals_total", "n_arrays", "frac_nonzero",
                       "frac_delta_gt0", "p75_delta", "p90_delta", "null_p95",
                       "frac_fg_above_null_p95"))
  fwrite(zmeas, file.path(RESULTS_DIR, paste0("period_enrichment_", OUT_SUFFIX, "_zeroinfl_measures.csv")))
  message("Saved: period_enrichment_merged_zeroinfl_measures.csv")

  zplot <- zmeas[sample %in% c("XG_150", "XG_151")]
  zplot[, sample_label := factor(sample_label, levels = c("CENP-A rep1", "CENP-A rep2"))]
  zlong <- melt(zplot,
    measure.vars = c("frac_delta_gt0", "frac_fg_above_null_p95"),
    variable.name = "measure", value.name = "frac")
  zlong[, measure := factor(measure,
    levels = c("frac_delta_gt0", "frac_fg_above_null_p95"),
    labels = c("Fraction of arrays above null baseline (Δ > 0)",
               "Fraction of arrays above null 95th percentile"))]

  p_zinfl <- ggplot(zlong,
                    aes(x = factor(bin_id, levels = TESTED_BINS), y = frac,
                        color = sample_label, group = sample_label)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.5) +
    facet_wrap(~ measure, ncol = 1, scales = "free_y") +
    scale_color_manual(values = SAMPLE_COLORS, name = NULL) +
    scale_x_discrete(labels = BIN_SHORT[as.character(TESTED_BINS)]) +
    labs(x = "TRF repeat period bin", y = "Fraction of observed arrays") +
    theme_period +
    theme(axis.text.x = element_text(size = FONT_SIZE - 3))

  ggsave(file.path(PLOTS_DIR, paste0("supp_zeroinflation_measures_", OUT_SUFFIX, ".pdf")),
         p_zinfl, width = 8, height = 6)
  ggsave(file.path(PLOTS_DIR, paste0("supp_zeroinflation_measures_", OUT_SUFFIX, ".png")),
         p_zinfl, width = 8, height = 6, dpi = 300)
  message("Saved: supp_zeroinflation_measures_merged.pdf/png")
}

# ============================================================================
# SECTION 8: Save results
# ============================================================================
message("\n=== SECTION 8: Saving results ===")

fwrite(signal, file.path(RESULTS_DIR, paste0("period_enrichment_", OUT_SUFFIX, "_foreground.csv")))

if (has_bg) {
  fwrite(perm_pvals, file.path(RESULTS_DIR, paste0("period_enrichment_", OUT_SUFFIX, "_permutation_results.csv")))

  sink(file.path(RESULTS_DIR, paste0("period_enrichment_", OUT_SUFFIX, "_permutation_results.txt")))
  cat("TRF Period-Size Stratified CENP-A Enrichment — Permutation Test Results (MERGED ARRAYS)\n")
  cat("========================================================================\n\n")
  cat(sprintf("Date: %s\n", Sys.time()))
  cat(sprintf("Unit of observation: MERGED REPEAT ARRAYS (bedtools merge -d 0 per period bin)\n"))
  cat(sprintf("Background: %d iterations per bin, chr- and length-matched shuffle\n",
              max(perm_pvals$n_iters)))
  cat(sprintf("Bins tested: %s\n", paste(TESTED_BINS, collapse = ", ")))
  cat(sprintf("P-values: empirical, MATCHED test. Bins 1-5 use paired 5,000-array\n"))
  cat(sprintf("          resampling: observed and shuffled medians are computed on the\n"))
  cat(sprintf("          SAME merged arrays per iteration, D = observed_median - shuffled_median.\n"))
  cat(sprintf("          Bins 6-9 shuffle the full array set (observed = full-bin median).\n"))
  cat(sprintf("          Enrichment P = fraction of iterations with D <= 0.\n"))
  cat(sprintf("FDR correction: Benjamini-Hochberg per sample\n\n"))

  for (s in c("XG_150", "XG_151", "XG_152", "XG_153")) {
    cat(sprintf("=== %s ===\n", SAMPLE_LABELS[s]))
    sub <- perm_pvals[sample == s][order(bin_id)]
    if (nrow(sub) == 0) { cat("  No data\n\n"); next }
    cat(sprintf("%-6s %-30s %8s %8s %8s %8s %8s %8s %8s %8s\n",
                "Bin", "Label", "Fg_med", "Bg_null", "Bg_SD", "D_mean", "D_SD", "Z", "P_enr", "P_adj"))
    cat(sprintf("%-6s %-30s %8s %8s %8s %8s %8s %8s %8s %8s\n",
                "---", "-----", "------", "-------", "-----", "-------", "-----", "-", "-----", "-----"))
    for (i in seq_len(nrow(sub))) {
      r <- sub[i]
      zs <- if (is.na(r$z_score)) "    —" else sprintf("%8.2f", r$z_score)
      cat(sprintf("%-6s %-30s %8.3f %8.3f %8.3f %8.3f %8.3f %8s %8.3f %8.3f %s\n",
                  r$bin_id, r$bin_label, r$fg_median,
                  r$bg_mean_of_medians, r$bg_sd_of_medians,
                  r$paired_effect, r$paired_effect_sd,
                  zs, r$p_enrich, r$p_adj, r$sig))
    }
    cat("\n")
  }

  cat("\nRedundancy metrics (merged arrays, -d 0):\n")
  if (file.exists(stats_file)) {
    for (i in seq_len(nrow(array_stats))) {
      r <- array_stats[i]
      cat(sprintf("  bin %s (%s): %s intervals -> %s arrays | avg %.2f | median %d | p75 %d | max %d intervals/array | %s%% of intervals in multi arrays\n",
                  r$bin_id, r$bin_label, format(r$n_intervals_total, big.mark=","),
                  format(r$n_arrays, big.mark=","), r$avg_intervals_per_array,
                  r$median_intervals_per_array, r$p75_intervals_per_array,
                  r$max_intervals_per_array, sprintf("%.1f", 100*r$frac_intervals_in_multi)))
    }
  }

  cat("\nKruskal-Wallis (across all bins 1-9, merged arrays):\n")
  for (s in c("XG_150", "XG_151")) {
    kw <- kw_results[[s]]$kw
    cat(sprintf("  %s: chi2=%.2f, df=%d, P=%.3e\n", SAMPLE_LABELS[s], kw$statistic, kw$parameter, kw$p.value))
  }

  sink()
  message("\nResults saved to: ", RESULTS_DIR)
}

writeLines(capture.output(sessionInfo()), file.path(RESULTS_DIR, paste0("session_info_", OUT_SUFFIX, ".txt")))

message("\n=== Analysis complete (merged arrays) ===")
message("Plots: ", PLOTS_DIR)
message("Results: ", RESULTS_DIR)
