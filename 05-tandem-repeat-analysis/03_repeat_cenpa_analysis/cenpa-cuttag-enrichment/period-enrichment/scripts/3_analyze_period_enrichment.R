#!/usr/bin/env Rscript
# ============================================================================
# 3_analyze_period_enrichment.R
# Period-size stratified CENP-A CUT&Tag enrichment at TRF tandem repeat loci.
#
# Primary inference: permutation test (foreground vs chromosome- and length-
# matched shuffled intervals, 1,000 iterations).
# Secondary: Kruskal-Wallis across all bins.
#
# Bin scheme (2026-08-03, 9 bins; monomer families 193-195/348-349/386-390 isolated):
#   1 = 1-10 bp (microsatellites)     6 = 348-349 bp
#   2 = 11-50 bp (minisatellites)     7 = 350-385 bp
#   3 = 51-192 bp                     8 = 386-390 bp
#   4 = 193-195 bp                    9 = 391+ bp
#   5 = 196-347 bp
# Bins 1-5 (many intervals) use a 5,000-interval subsample null; bins 6-9 (tiny) full-set.
#
# Usage:
#   conda activate r-visualizations
#   Rscript 3_analyze_period_enrichment.R
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
COUNTS_DIR  <- file.path(DATA_DIR, "counts")
BG_COUNTS_DIR <- file.path(DATA_DIR, "permutation", "counts")
PLOTS_DIR   <- file.path(PERIOD_DIR, "plots")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

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

# Period bin definitions (2026-08-03, 9 bins; monomer families isolated)
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

# Bins tested by the permutation test (all 9 analyzed bins)
TESTED_BINS <- 1:9
# Full-set null bins (tiny interval counts; no subsampling)
CENTROMERIC_BINS <- c(6, 7, 8, 9)
# Subsample-null bins (many intervals -> 5,000-interval subsample)
CONTROL_BINS <- c(1, 2, 3, 4, 5)

# ============================================================================
# ggplot2 theme
# ============================================================================
theme_period <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

# ============================================================================
# SECTION 1: Load foreground data
# ============================================================================
message("=== SECTION 1: Loading foreground data ===")

# Library sizes
lib_file <- file.path(PARENT_DIR, "data", "qc", "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
} else {
  frag_dir <- file.path(PARENT_DIR, "data", "fragments")
  lib_sizes <- sapply(SAMPLES, function(s) {
    f <- file.path(frag_dir, paste0(s, "_fragments.bed"))
    if (file.exists(f)) as.numeric(system(paste("wc -l <", f), intern = TRUE)) else NA
  })
  names(lib_sizes) <- SAMPLES
}
message("Library sizes (fragments):")
for (s in SAMPLES) message(sprintf("  %s: %d", s, lib_sizes[s]))

# Load foreground signal
message("\nLoading foreground signal...")
signal_list <- list()
for (s in SAMPLES) {
  f <- file.path(COUNTS_DIR, paste0("trf_signal_", s, ".tsv"))
  if (!file.exists(f)) {
    message(sprintf("  WARNING: Signal file not found: %s", f))
    next
  }
  x <- fread(f, header = FALSE,
             col.names = c("chrom", "start", "end", "interval_id",
                           "bin_id", "period_size", "copies_aligned",
                           "match_percent", "mean_coverage"))
  x[, sample := s]
  x[, interval_length := end - start]
  signal_list[[s]] <- x
  message(sprintf("  %s: %d intervals", s, nrow(x)))
}
signal <- rbindlist(signal_list)
message(sprintf("  Total foreground: %d rows", nrow(signal)))

if (nrow(signal) == 0) {
  stop("No foreground signal data. Run 2_extract_signal.sh first.")
}

signal[, sample_label := factor(SAMPLE_LABELS[sample],
                                 levels = c("CENP-A rep1", "CENP-A rep2",
                                            "H3K27ac", "H3K27ac rep2"))]
signal[, bin_factor := factor(bin_id, levels = 1:9)]

# ============================================================================
# SECTION 2: Normalize
# ============================================================================
message("\n=== SECTION 2: Normalizing signal ===")

signal[, lib_size := lib_sizes[sample]]
signal[, norm_signal := mean_coverage / (lib_size / 1e6)]

pseudo <- min(signal$norm_signal[signal$norm_signal > 0], na.rm = TRUE) / 2
message(sprintf("Pseudocount: %.4e", pseudo))

signal[, log2_signal := log2(norm_signal + pseudo)]

# ============================================================================
# SECTION 3: Load background (permuted) data
# ============================================================================
message("\n=== SECTION 3: Loading background (permuted) data ===")

bg_list <- list()
has_bg <- FALSE
for (s in SAMPLES) {
  f <- file.path(BG_COUNTS_DIR, paste0("trf_bg_signal_", s, ".tsv"))
  if (!file.exists(f)) {
    message(sprintf("  WARNING: No background data for %s (run 2b_period_background.sh)", s))
    next
  }
  # Columns from bedtools map + our added columns:
  # chr, start, end, bin_id, iter, mean_coverage  (6 cols, no header)
  x <- fread(f, header = FALSE,
             col.names = c("chrom", "start", "end", "bin_id", "iter", "mean_coverage"))
  x[, sample := s]
  # Normalize per sample, then keep only the columns needed downstream so the
  # combined bg table (~111M rows) stays small enough to hold on one node.
  x[, norm_signal := mean_coverage / (lib_sizes[s] / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
  n_bg <- nrow(x)
  n_bins <- uniqueN(x$bin_id)
  bg_list[[s]] <- x[, .(sample, chrom, bin_id, iter, log2_signal)]
  has_bg <- TRUE
  message(sprintf("  %s: %d background intervals across %d bins", s, n_bg, n_bins))
}

# Load control-permutation nulls for bins 1-5 (micro/mini subsample).
# Different format: 7 cols WITH header (bin_id, interval_id, iter, chrom, start, end, mean_signal).
# The interval_id column lets us match the foreground median to the same 5,000-interval
# subsample that defines the null (required for a balanced fg-vs-null comparison).
ctrl_has_data <- FALSE
for (s in SAMPLES) {
  f <- file.path(COUNTS_DIR, paste0("trf_ctrl_bg_signal_", s, ".tsv"))
  if (!file.exists(f)) {
    message(sprintf("  NOTE: No control-permutation background for %s (run 2c_control_permutation.sh)", s))
    next
  }
  x <- fread(f, header = TRUE)  # bin_id, interval_id, iter, chrom, start, end, mean_signal
  setnames(x, "mean_signal", "mean_coverage")
  x[, sample := s]
  x[, interval_length := end - start]
  # Join OBSERVED signal for the same interval_ids, so the matched permutation
  # test can compare observed vs shuffled medians computed on the SAME 5,000
  # intervals in each iteration (the shuffle output preserves interval_id).
  obs <- signal[sample == s, .(interval_id, obs_coverage = mean_coverage)]
  x <- merge(x, obs, by = "interval_id", all.x = TRUE)
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
  message(sprintf("  %s: %d control-permutation background intervals across %d bins",
                  s, nrow(x), uniqueN(x$bin_id)))
}

if (has_bg) {
  bg <- rbindlist(bg_list, fill = TRUE)
  message(sprintf("  Total background: %d rows", nrow(bg)))
  message("  (compact: sample, chrom, bin_id, iter, log2_signal, log2_obs)")

  # Check iterations per bin
  bg_iters <- bg[, .(n_iters = uniqueN(iter)), by = .(sample, bin_id)]
  message("\nBackground iterations per bin per sample:")
  message(sprintf("  Range: %d - %d iterations", min(bg_iters$n_iters), max(bg_iters$n_iters)))
}

# ============================================================================
# SECTION 4: Permutation test (primary inference)
# ============================================================================
message("\n=== SECTION 4: Permutation test ===")

if (has_bg) {
  # Descriptive foreground summary (full-bin; kept for display/reporting).
  # All intervals in the bin — no subsampling on the descriptive side.
  fg_bin <- signal[bin_id %in% TESTED_BINS,
                    .(fg_median = median(log2_signal, na.rm = TRUE),
                      fg_mean = mean(log2_signal, na.rm = TRUE),
                      n_fg = .N),
                    by = .(sample, bin_id)]

  # Matched (paired) permutation: per-iteration observed AND shuffled medians
  # computed on the SAME interval set, so the observed statistic and each null
  # statistic use identical sample sizes and the same procedure.
  #   bins 1-5: obs_median = median observed log2 signal over the 5,000
  #             intervals that were subsampled AND shuffled in that iteration
  #             (ctrl null rows carry interval_id; observed coverage joined in
  #             Section 3). The null marginalizes over fresh 5,000-subsample
  #             draws per iteration (2d_permutation_per_bin.sh).
  #   bins 6-9: shuffled set = full interval set, so the observed side is the
  #             full-bin median (constant across iterations).
  bg_bin_iter <- bg[bin_id %in% TESTED_BINS,
                     .(bg_median = median(log2_signal, na.rm = TRUE),
                       bg_mean = mean(log2_signal, na.rm = TRUE),
                       obs_median = median(log2_obs, na.rm = TRUE)),
                     by = .(sample, bin_id, iter)]
  bg_bin_iter <- merge(bg_bin_iter, fg_bin[, .(sample, bin_id, fg_median)],
                       by = c("sample", "bin_id"))
  bg_bin_iter[, obs_med := ifelse(is.na(obs_median), fg_median, obs_median)]

  # Paired matched effect per iteration: D_j = observed median - shuffled median
  bg_bin_iter[, D := obs_med - bg_median]

  # Empirical P-value from the matched D_j distribution.
  #   p_enrich  = fraction of iterations where observed <= shuffled (D <= 0)
  #   p_deplete = fraction of iterations where observed >= shuffled (D >= 0)
  perm_pvals <- bg_bin_iter[, .(
    fg_median = fg_median[1],
    p_enrich = mean(D <= 0),
    p_deplete = mean(D >= 0),
    n_iters = .N,
    bg_mean_of_medians = mean(bg_median),
    bg_sd_of_medians = sd(bg_median),
    obs_mean_of_medians = mean(obs_med),
    paired_effect = mean(D),        # mean matched observed-minus-null effect
    paired_effect_sd = sd(D),
    z_score = mean(D) / (sd(D) + 1e-10)
  ), by = .(sample, bin_id)]
  # Mask degenerate z-scores. When the null is a point mass (every iteration
  # yields the same median — e.g. zero-inflated bins 1-3 where >50% of shuffled
  # positions have no coverage), the matched effect D is also constant, SD = 0
  # and z = mean(D)/0 is meaningless. The empirical p-values remain valid; only
  # the z column is masked.
  perm_pvals[paired_effect_sd == 0, z_score := NA_real_]

  # Add bin labels
  perm_pvals[, bin_label := BIN_LABELS_FULL[as.character(bin_id)]]
  perm_pvals[, sample_label := SAMPLE_LABELS[sample]]

  # Two-sided P-value
  perm_pvals[, p_two_sided := pmin(p_enrich, p_deplete) * 2]
  perm_pvals[, p_two_sided := pmin(p_two_sided, 1.0)]

  # Apply FDR correction per sample
  perm_pvals[, p_adj := p.adjust(p_two_sided, method = "BH"), by = sample]

  # Significance annotation
  perm_pvals[, sig := ifelse(p_adj < 0.001, "***",
                      ifelse(p_adj < 0.01, "**",
                      ifelse(p_adj < 0.05, "*",
                      ifelse(p_adj < 0.1, ".", "ns"))))]

  message(sprintf("\nperm_pvals dimensions: %d rows x %d cols", nrow(perm_pvals), ncol(perm_pvals)))
  message(sprintf("perm_pvals columns: %s", paste(names(perm_pvals), collapse = ", ")))
  message(sprintf("perm_pvals samples: %s", paste(unique(perm_pvals$sample), collapse = ", ")))

  message("\nPermutation test results (empirical P, one-sided enrichment):")
  message("===============================================================")
  for (s in c("XG_150", "XG_151")) {
    message(sprintf("\n%s:", SAMPLE_LABELS[s]))
    sub <- perm_pvals[sample == s][order(bin_id)]
    message(sprintf("  (sub-table has %d rows)", nrow(sub)))
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
  message("  SKIP: No background data available. Run 2b_period_background.sh first.")
}

# ============================================================================
# SECTION 5: Kruskal-Wallis (secondary — across-bin comparison)
# ============================================================================
message("\n=== SECTION 5: Kruskal-Wallis ===")

kw_results <- list()
for (s in c("XG_150", "XG_151")) {
  sub <- signal[sample == s & bin_id <= 9]
  kw <- kruskal.test(log2_signal ~ bin_factor, data = sub)
  message(sprintf("%s: chi2=%.2f, df=%d, P=%.3e",
                  SAMPLE_LABELS[s], kw$statistic, kw$parameter, kw$p.value))

  # Dunn's post-hoc
  dunn_res <- sub %>%
    dunn_test(log2_signal ~ bin_factor, p.adjust.method = "bonferroni") %>%
    as.data.table()
  kw_results[[s]] <- list(kw = kw, dunn = dunn_res)
}

# ============================================================================
# SECTION 6: Main figure — Period-stratified enrichment with permutation null
# ============================================================================
message("\n=== SECTION 6: Main figure ===")

# Panel A: Foreground vs background per bin (CENP-A samples, bins 6-9)
if (has_bg) {
  # Prepare foreground bin summary for bins 6-9
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

  # Prepare background distribution for bins 6-9
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
    # Background null distribution (grey)
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
    # Foreground median (colored points)
    geom_point(data = fg_bin_plot,
               aes(x = as.numeric(as.character(bin_id)), y = fg_median,
                   color = sample_label),
               size = 3, shape = 18) +
    # Significance labels
    geom_text(data = fg_bin_plot[sig != " "],
              aes(x = as.numeric(as.character(bin_id)),
                  y = fg_median + 0.2, label = sig),
              size = 5, color = "red", fontface = "bold") +
    # Background null median reference
    geom_point(data = bg_bin_dist,
               aes(x = as.numeric(as.character(bin_id)), y = median_of_median),
               shape = 3, size = 1.5, color = "grey40") +
    facet_wrap(~ sample_label, ncol = 1) +
    scale_color_manual(values = SAMPLE_COLORS, guide = "none") +
    scale_x_continuous(
      breaks = TESTED_BINS,
      labels = BIN_SHORT[as.character(TESTED_BINS)]
    ) +
    labs(
      x = "TRF repeat period bin",
      y = expression("Median log"[2] * "(normalized CUT&Tag signal)"),
      title = "CENP-A enrichment at TRF repeat loci — permutation test",
      subtitle = paste0("Grey: 95% and 50% CI of null distribution (1,000 shuffles). ",
                        "Filled diamond: observed foreground median. ",
                        "Cross: null median.\n",
                        "Bins 1-5: matched paired test on 5,000-interval subsamples; ",
                        "bins 6-9: full set. *** P<0.001, ** P<0.01, * P<0.05, . P<0.1 (FDR-corrected)")
    ) +
    theme_period

  p_perm_full <- p_perm + theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5))
} else {
  # Fallback: boxplot as before
  message("No background data — using boxplot-only figure")
  cenpa_signal <- signal[sample %in% c("XG_150", "XG_151") & bin_id <= 9]

  p_perm_full <- ggplot(cenpa_signal, aes(x = bin_factor, y = log2_signal)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_boxplot(aes(fill = bin_id %in% TESTED_BINS),
                 outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.25) +
    facet_wrap(~ sample_label, ncol = 1) +
    scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "grey80"), guide = "none") +
    scale_x_discrete(labels = BIN_SHORT[as.character(1:10)]) +
    labs(x = "TRF repeat period", y = expression(log[2] * "(normalized signal)"),
         title = "CENP-A enrichment by TRF repeat period",
         subtitle = "Boxplot of per-interval signal (fallback without permutation null).") +
    theme_period +
    theme(axis.text.x = element_text(size = 7))
}

ggsave(file.path(PLOTS_DIR, "period_permutation_enrichment.pdf"),
       p_perm_full, width = 8, height = 7)
ggsave(file.path(PLOTS_DIR, "period_permutation_enrichment.png"),
       p_perm_full, width = 8, height = 7, dpi = 300)
message("Saved: period_permutation_enrichment.pdf/png")

# Panel B: All bins overview (foreground only, CENP-A + H3K27ac comparison)
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
  scale_x_continuous(breaks = 1:10, labels = BIN_SHORT[as.character(1:10)]) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(x = "TRF repeat period", y = expression("Mean log"[2] * "(normalized signal) ± SE"),
       title = "Mean CENP-A vs H3K27ac signal across all period bins",
       subtitle = "Error bars: ±1 SE. Point size ∝ number of intervals.") +
  theme_period +
  theme(axis.text.x = element_text(size = 7))

ggsave(file.path(PLOTS_DIR, "period_all_bins_overview.pdf"),
       p_overview, width = 10, height = 5)
ggsave(file.path(PLOTS_DIR, "period_all_bins_overview.png"),
       p_overview, width = 10, height = 5, dpi = 300)
message("Saved: period_all_bins_overview.pdf/png")

# ============================================================================
# SECTION 7: Supplementary figures
# ============================================================================
message("\n=== SECTION 7: Supplementary figures ===")

# Supp A: Z-score heatmap per chromosome per bin (CENP-A rep1)
if (has_bg) {
  # Compute per-chromosome per-bin foreground median (ALL intervals in each bin)
  chr_fg <- signal[sample == "XG_150" & bin_id %in% TESTED_BINS,
                    .(fg_median = median(log2_signal, na.rm = TRUE)),
                    by = .(chrom, bin_id)]

  # Compute per-chromosome per-bin background null
  chr_bg <- bg[sample == "XG_150" & bin_id %in% TESTED_BINS,
                .(bg_median = median(log2_signal, na.rm = TRUE)),
                by = .(chrom, bin_id, iter)]
  chr_bg_null <- chr_bg[, .(
    null_mean = mean(bg_median),
    null_sd = sd(bg_median)
  ), by = .(chrom, bin_id)]

  # Z-scores (mask degenerate cells where the per-chromosome null is a point mass, null_sd = 0)
  chr_z <- merge(chr_fg, chr_bg_null, by = c("chrom", "bin_id"))
  chr_z[, z_score := (fg_median - null_mean) / (null_sd + 1e-10)]
  chr_z[null_sd == 0, z_score := NA_real_]
  chr_z[, z_capped := ifelse(is.na(z_score), NA_real_, pmax(pmin(z_score, 5), -5))]

  # Filter chromosomes with data in bin 7
  chr_with_data <- chr_z[bin_id == 7, unique(chrom)]
  chr_z <- chr_z[chrom %in% chr_with_data]
  chr_order <- chr_z[bin_id == 7][order(-z_score)]$chrom
  chr_z[, chrom_factor := factor(chrom, levels = chr_order)]
  chr_z[, bin_factor := factor(bin_id, levels = TESTED_BINS)]

  p_zheat <- ggplot(chr_z, aes(x = bin_factor, y = chrom_factor, fill = z_capped)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is.na(z_score), "", sprintf("%.1f", z_score))),
              size = 2.8) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                          midpoint = 0, na.value = "grey90", name = "Z-score") +
    scale_x_discrete(labels = BIN_SHORT[as.character(TESTED_BINS)], position = "top") +
    labs(x = "TRF period bin", y = "Chromosome",
         title = "CENP-A enrichment Z-scores by bin and chromosome (rep1)",
         subtitle = "Z = (observed median - null mean) / null SD. Red = enriched, Blue = depleted.") +
    theme_period +
    theme(panel.grid = element_blank(), legend.position = "right")

  ggsave(file.path(PLOTS_DIR, "supp_permutation_zscore_heatmap.pdf"),
         p_zheat, width = 9, height = 8)
  ggsave(file.path(PLOTS_DIR, "supp_permutation_zscore_heatmap.png"),
         p_zheat, width = 9, height = 8, dpi = 300)
  message("Saved: supp_permutation_zscore_heatmap.pdf/png")
}

# Supp B: 348-349 bp monomer (bin 6) null distribution detail.
# NOTE: in the 9-bin scheme the 349-bp monomer family is bin 6 (348-349 bp),
# not bin 7 (350-385 bp). Was `bin_id == 7` before the 9-bin rebin (2026-08-02).
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
         y = "Frequency (1,000 iterations)",
         title = "348-349 bp monomer (bin 6): Null distribution vs observed CENP-A",
         subtitle = "Histogram: null medians. Vertical line: observed foreground median.") +
    theme_period

  ggsave(file.path(PLOTS_DIR, "supp_bin6_null_distribution.pdf"),
         p_null, width = 7, height = 6)
  ggsave(file.path(PLOTS_DIR, "supp_bin6_null_distribution.png"),
         p_null, width = 7, height = 6, dpi = 300)
  message("Saved: supp_bin6_null_distribution.pdf/png")
} else {
  message("  SKIP supp_bin6: perm_pvals missing or empty")
}

# Supp C: Copy number vs CENP-A signal (within 348-349 bp monomer, bin 6)
bin6_signal <- signal[bin_id == 6 & sample %in% c("XG_150", "XG_151")]
bin6_signal[, sample_label := factor(SAMPLE_LABELS[sample],
                                      levels = c("CENP-A rep1", "CENP-A rep2"))]

p_copies <- ggplot(bin6_signal, aes(x = log10(copies_aligned + 1), y = log2_signal)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_point(aes(color = chrom), size = 0.5, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "#2166AC", linewidth = 0.8) +
  facet_wrap(~ sample_label, nrow = 1) +
  scale_color_discrete(guide = "none") +
  labs(x = expression(log[10] * "(copy number + 1)"),
       y = expression(log[2] * "(normalized signal)"),
       title = "348-349 bp monomer: copy number vs CENP-A signal") +
  theme_period

ggsave(file.path(PLOTS_DIR, "supp_349bp_copies_vs_signal.pdf"), p_copies, width = 9, height = 4)
ggsave(file.path(PLOTS_DIR, "supp_349bp_copies_vs_signal.png"), p_copies, width = 9, height = 4, dpi = 300)
message("Saved: supp_349bp_copies_vs_signal.pdf/png")

# Supp D: Match% (repeat homogeneity) vs CENP-A signal
p_match <- ggplot(bin6_signal, aes(x = match_percent, y = log2_signal)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_point(aes(color = chrom), size = 0.5, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "#2166AC", linewidth = 0.8) +
  facet_wrap(~ sample_label, nrow = 1) +
  scale_color_discrete(guide = "none") +
  labs(x = "Match percent (%)", y = expression(log[2] * "(normalized signal)"),
       title = "348-349 bp monomer: repeat homogeneity vs CENP-A signal") +
  theme_period

ggsave(file.path(PLOTS_DIR, "supp_349bp_match_vs_signal.pdf"), p_match, width = 9, height = 4)
ggsave(file.path(PLOTS_DIR, "supp_349bp_match_vs_signal.png"), p_match, width = 9, height = 4, dpi = 300)
message("Saved: supp_349bp_match_vs_signal.pdf/png")

# ============================================================================
# SECTION 7b: Zero-inflation & non-independence supplement
# ============================================================================
message("\n=== SECTION 7b: Zero-inflation & non-independence supplement ===")

# ── Repeat-array counts (non-independence): merge overlapping intervals per bin ──
TRF_BED <- file.path(DATA_DIR, "trf_chr_only_period_bins.bed")
if (file.exists(TRF_BED)) {
  trf <- fread(TRF_BED, header = FALSE,
    col.names = c("chrom", "start", "end", "interval_id", "bin_id",
                  "period", "copies", "match_pct"))
  arr <- trf[bin_id %in% TESTED_BINS][order(chrom, start, end)]
  arr[, prev_end := shift(end, type = "lag"), by = .(chrom, bin_id)]
  arr[, new_array := is.na(prev_end) | start > prev_end]
  arr[, array_id := cumsum(new_array), by = .(chrom, bin_id)]
  n_arrays <- arr[, .(n_arrays = uniqueN(array_id)), by = bin_id]
  message("Repeat-array counts per bin (merged overlapping intervals):")
  for (b in TESTED_BINS) {
    message(sprintf("  bin %d: %d arrays", b, n_arrays[bin_id == b, n_arrays]))
  }
} else {
  n_arrays <- data.table(bin_id = TESTED_BINS, n_arrays = NA_integer_)
  message("  WARNING: TRF bed not found — n_arrays set to NA")
}

if (has_bg) {
  # Null baseline (mean of per-iteration shuffled medians; same as the violin)
  null_base <- bg_bin_iter[, .(null_median = mean(bg_median)), by = .(sample, bin_id)]

  fg_all <- signal[bin_id %in% TESTED_BINS]
  fg_all <- merge(fg_all, null_base, by = c("sample", "bin_id"))
  fg_all[, delta := log2_signal - null_median]

  # Pooled per-interval null 95th percentile (bins 1-5: shuffled ctrl intervals;
  # bins 6-9: full-set shuffled intervals)
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
  zmeas <- merge(zmeas, n_arrays, by = "bin_id")
  zmeas <- merge(zmeas, null_p95, by = c("sample", "bin_id"))
  zmeas[, bin_label := BIN_LABELS_FULL[as.character(bin_id)]]
  zmeas[, sample_label := SAMPLE_LABELS[sample]]
  setcolorder(zmeas, c("sample", "sample_label", "bin_id", "bin_label", "n",
                       "n_arrays", "frac_nonzero", "frac_delta_gt0",
                       "p75_delta", "p90_delta", "null_p95",
                       "frac_fg_above_null_p95"))
  fwrite(zmeas, file.path(RESULTS_DIR, "period_enrichment_zeroinfl_measures.csv"))
  message("Saved: period_enrichment_zeroinfl_measures.csv")

  # Compact supplementary figure (CENP-A replicates): subset-enrichment measures
  zplot <- zmeas[sample %in% c("XG_150", "XG_151")]
  zplot[, sample_label := factor(sample_label,
                                  levels = c("CENP-A rep1", "CENP-A rep2"))]
  zlong <- melt(zplot,
    measure.vars = c("frac_delta_gt0", "frac_fg_above_null_p95"),
    variable.name = "measure", value.name = "frac")
  zlong[, measure := factor(measure,
    levels = c("frac_delta_gt0", "frac_fg_above_null_p95"),
    labels = c("Fraction of intervals above null baseline (Δ > 0)",
               "Fraction of intervals above null 95th percentile"))]

  p_zinfl <- ggplot(zlong,
                    aes(x = factor(bin_id, levels = TESTED_BINS), y = frac,
                        color = sample_label, group = sample_label)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.5) +
    facet_wrap(~ measure, ncol = 1, scales = "free_y") +
    scale_color_manual(values = SAMPLE_COLORS, name = NULL) +
    scale_x_discrete(labels = BIN_SHORT[as.character(TESTED_BINS)]) +
    labs(x = "TRF repeat period bin", y = "Fraction of observed intervals",
         title = "Zero-inflated TRF bins: subset-enrichment measures (CENP-A)",
         subtitle = "Bins 1-3 sit near the pseudocount floor; these measures distinguish a minority of enriched loci from uniformly inactive loci.") +
    theme_period +
    theme(axis.text.x = element_text(size = 7))

  ggsave(file.path(PLOTS_DIR, "supp_zeroinflation_measures.pdf"), p_zinfl,
         width = 8, height = 6)
  ggsave(file.path(PLOTS_DIR, "supp_zeroinflation_measures.png"), p_zinfl,
         width = 8, height = 6, dpi = 300)
  message("Saved: supp_zeroinflation_measures.pdf/png")
}

# ============================================================================
# SECTION 8: Save results
# ============================================================================
message("\n=== SECTION 8: Saving results ===")

# Full per-interval foreground data
fwrite(signal, file.path(RESULTS_DIR, "period_enrichment_foreground.csv"))

# Permutation results
if (has_bg) {
  fwrite(perm_pvals, file.path(RESULTS_DIR, "period_enrichment_permutation_results.csv"))

  sink(file.path(RESULTS_DIR, "period_enrichment_permutation_results.txt"))
  cat("TRF Period-Size Stratified CENP-A Enrichment — Permutation Test Results\n")
  cat("========================================================================\n\n")
  cat(sprintf("Date: %s\n", Sys.time()))
  cat(sprintf("Background: %d iterations per bin, chr- and length-matched shuffle\n",
              max(perm_pvals$n_iters)))
  cat(sprintf("Bins tested: %s\n", paste(TESTED_BINS, collapse = ", ")))
  cat(sprintf("P-values: empirical, MATCHED test. Bins 1-5 use paired 5,000-interval\n"))
  cat(sprintf("          resampling: observed and shuffled medians are computed on the\n"))
  cat(sprintf("          SAME intervals per iteration, D = observed_median - shuffled_median.\n"))
  cat(sprintf("          Bins 6-9 shuffle the full interval set (observed = full-bin median).\n"))
  cat(sprintf("          Enrichment P = fraction of iterations with D <= 0.\n"))
  cat(sprintf("FDR correction: Benjamini-Hochberg per sample\n\n"))

  for (s in c("XG_150", "XG_151", "XG_152", "XG_153")) {
    cat(sprintf("=== %s ===\n", SAMPLE_LABELS[s]))
    sub <- perm_pvals[sample == s][order(bin_id)]
    if (nrow(sub) == 0) {
      cat("  No data\n\n")
      next
    }
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

  cat("\nKruskal-Wallis (across all bins 1-10):\n")
  for (s in c("XG_150", "XG_151")) {
    kw <- kw_results[[s]]$kw
    cat(sprintf("  %s: chi2=%.2f, df=%d, P=%.3e\n",
                SAMPLE_LABELS[s], kw$statistic, kw$parameter, kw$p.value))
  }

  sink()
  message("\nResults saved to: ", RESULTS_DIR)
}

# Session info
writeLines(capture.output(sessionInfo()), file.path(RESULTS_DIR, "session_info.txt"))

message("\n=== Analysis complete ===")
message("Plots: ", PLOTS_DIR)
message("Results: ", RESULTS_DIR)
