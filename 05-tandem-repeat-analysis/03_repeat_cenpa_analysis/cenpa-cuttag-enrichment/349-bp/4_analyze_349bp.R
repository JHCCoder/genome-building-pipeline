#!/usr/bin/env Rscript
# ============================================================================
# 4_analyze_349bp.R — CENP-A enrichment at 349-bp tandem repeat / HOR loci
#
# Tests whether ~349-bp repeat loci (defined independently by TRF + centroAnno
# HORs) carry more CENP-A CUT&Tag signal than chromosome- and length-matched
# randomized intervals.
#
# Primary inference: per-chromosome median log2(normalized signal) delta (fg - bg)
#   Paired Wilcoxon signed-rank test across chromosomes (H0: delta = 0)
#
# Usage:
#   conda activate r-visualizations
#   Rscript 4_analyze_349bp.R
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
BASE_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/349-bp"
PARENT_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
DATA_DIR <- file.path(BASE_DIR, "data")
COUNTS_DIR <- file.path(DATA_DIR, "counts")
PLOTS_DIR <- file.path(BASE_DIR, "plots")
RESULTS_DIR <- file.path(BASE_DIR, "results")

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Chromosomes — chrY excluded (no 349-bp TRF repeats)
CHROMOSOMES <- paste0("chr", c(1:28, "X"))

SAMPLES <- c("XG_150", "XG_151", "XG_152", "XG_153")
SAMPLE_LABELS <- c(
  "XG_150" = "CENP-A rep1",
  "XG_151" = "CENP-A rep2",
  "XG_152" = "H3K27ac",
  "XG_153" = "H3K27ac rep2"
)
SAMPLE_COLORS <- c(
  "XG_150" = "#2166AC",
  "XG_151" = "#92C5DE",
  "XG_152" = "#B2182B",
  "XG_153" = "#D6604D"
)

# ============================================================================
# ggplot2 theme
# ============================================================================
theme_349 <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

# ============================================================================
# SECTION 1: Load data
# ============================================================================
message("=== SECTION 1: Loading data ===")

# Library sizes (from parent pipeline)
lib_file <- file.path(PARENT_DIR, "data", "qc", "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
} else {
  # Compute from fragment BEDs
  message("Computing library sizes from fragment BEDs...")
  frag_dir <- file.path(PARENT_DIR, "data", "fragments")
  lib_sizes <- sapply(SAMPLES, function(s) {
    f <- file.path(frag_dir, paste0(s, "_fragments.bed"))
    if (file.exists(f)) as.numeric(system(paste("wc -l <", f), intern = TRUE)) else NA
  })
  names(lib_sizes) <- SAMPLES
}
message("Library sizes (fragments):")
for (s in SAMPLES) message(sprintf("  %s: %d", s, lib_sizes[s]))

# Load foreground counts
message("\nLoading foreground counts...")
load_fg <- function(sample) {
  f <- file.path(COUNTS_DIR, paste0(sample, "_349bp_foreground.txt"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, header = FALSE, col.names = c("chrom", "start", "end", "interval_id", "size", "count"))
  x[, sample := sample]
  return(x)
}
fg_list <- lapply(SAMPLES, load_fg)
fg <- rbindlist(fg_list[!sapply(fg_list, is.null)])
message("  Foreground intervals: ", uniqueN(fg$interval_id), " per sample")

# Load background counts
message("Loading background counts...")
load_bg <- function(sample) {
  f <- file.path(COUNTS_DIR, paste0(sample, "_349bp_bg_shuffled.txt"))
  if (!file.exists(f)) return(NULL)
  x <- fread(f, header = FALSE,
             col.names = c("chrom", "start", "end", "iter_id", "interval_id", "count"))
  x[, sample := sample]
  return(x)
}
bg_list <- lapply(SAMPLES, load_bg)
bg <- rbindlist(bg_list[!sapply(bg_list, is.null)])
n_iter <- uniqueN(bg$iter_id)
message("  Background intervals: ", nrow(bg) / length(SAMPLES), " per sample (", n_iter, " iterations)")

# ============================================================================
# SECTION 2: Normalize CUT&Tag signal
# ============================================================================
message("\n=== SECTION 2: Normalize CUT&Tag signal ===")

# Add interval sizes to background (bg intervals are length-matched to fg intervals)
fg_sizes <- unique(fg[, .(interval_id, size)])
bg <- merge(bg, fg_sizes, by = "interval_id", all.x = TRUE)

compute_norm_signal <- function(dt, lib_sizes) {
  # dt must have: sample, size (interval length), count
  dt[, lib_size := lib_sizes[sample]]
  dt[, norm_signal := (count / (size / 1000)) / (lib_size / 1e6)]
  return(dt)
}

fg <- compute_norm_signal(fg, lib_sizes)
bg <- compute_norm_signal(bg, lib_sizes)

# Pseudocount: half the minimum non-zero normalized signal
pseudo_fg <- min(fg$norm_signal[fg$norm_signal > 0], na.rm = TRUE) / 2
pseudo_bg <- min(bg$norm_signal[bg$norm_signal > 0], na.rm = TRUE) / 2
pseudocount <- min(pseudo_fg, pseudo_bg)
message(sprintf("Pseudocount: %.2e", pseudocount))

fg[, log2_signal := log2(norm_signal + pseudocount)]
bg[, log2_signal := log2(norm_signal + pseudocount)]

# Add sample labels
fg[, sample_label := SAMPLE_LABELS[sample]]
bg[, sample_label := SAMPLE_LABELS[sample]]

message("Foreground: ", nrow(fg), " intervals across ", length(SAMPLES), " samples")
message("Background: ", nrow(bg), " intervals")

# ============================================================================
# SECTION 3: Primary inference — chromosome-level enrichment
# ============================================================================
message("\n=== SECTION 3: Primary inference — Chromosome-level enrichment ===")

# Per-chromosome foreground medians
chr_fg <- fg[, .(
  median_fg = median(log2_signal, na.rm = TRUE),
  mean_fg = mean(log2_signal, na.rm = TRUE),
  n_intervals = .N
), by = .(sample, chrom)]

# Per-chromosome background medians (aggregate across all iterations)
chr_bg <- bg[, .(
  median_bg = median(log2_signal, na.rm = TRUE),
  mean_bg = mean(log2_signal, na.rm = TRUE)
), by = .(sample, chrom)]

# Merge and compute delta
chr_delta <- merge(chr_fg, chr_bg, by = c("sample", "chrom"))
chr_delta[, delta := median_fg - median_bg]
chr_delta[, chrom_factor := factor(chrom, levels = CHROMOSOMES)]

# Chromosome ordering by mean CENPA delta
chr_order_data <- chr_delta[sample %in% c("XG_150", "XG_151"),
                             .(mean_delta = mean(delta, na.rm = TRUE)), by = chrom]
setorder(chr_order_data, mean_delta)
chr_delta[, chrom_ordered := factor(chrom, levels = chr_order_data$chrom)]
chr_delta[, sample_label := SAMPLE_LABELS[sample]]

# Statistical tests (Wilcoxon signed-rank)
chr_delta_wide <- dcast(chr_delta[sample %in% c("XG_150", "XG_151")],
                        chrom ~ sample, value.var = "delta")

wilcox_150 <- tryCatch(
  wilcox.test(chr_delta_wide$XG_150, mu = 0, alternative = "two.sided", exact = TRUE),
  error = function(e) wilcox.test(chr_delta_wide$XG_150, mu = 0, alternative = "two.sided", exact = FALSE)
)
wilcox_151 <- tryCatch(
  wilcox.test(chr_delta_wide$XG_151, mu = 0, alternative = "two.sided", exact = TRUE),
  error = function(e) wilcox.test(chr_delta_wide$XG_151, mu = 0, alternative = "two.sided", exact = FALSE)
)

message("Primary inference: Paired Wilcoxon signed-rank across chromosomes")
message(sprintf("  XG_150: V = %.0f, P = %.3e", wilcox_150$statistic, wilcox_150$p.value))
message(sprintf("  XG_151: V = %.0f, P = %.3e", wilcox_151$statistic, wilcox_151$p.value))

prop_pos_150 <- mean(chr_delta_wide$XG_150 > 0, na.rm = TRUE)
prop_pos_151 <- mean(chr_delta_wide$XG_151 > 0, na.rm = TRUE)
n_chr <- sum(!is.na(chr_delta_wide$XG_150))
message(sprintf("  Positive chromosomes: XG_150 = %d/%d (%.1f%%)",
                sum(chr_delta_wide$XG_150 > 0, na.rm = TRUE), n_chr, prop_pos_150 * 100))
message(sprintf("  Positive chromosomes: XG_151 = %d/%d (%.1f%%)",
                sum(chr_delta_wide$XG_151 > 0, na.rm = TRUE), n_chr, prop_pos_151 * 100))
message(sprintf("  Median delta (XG_150): %.3f", median(chr_delta_wide$XG_150, na.rm = TRUE)))
message(sprintf("  Median delta (XG_151): %.3f", median(chr_delta_wide$XG_151, na.rm = TRUE)))

# ============================================================================
# SECTION 4: Effect sizes
# ============================================================================
message("\n=== SECTION 4: Effect sizes ===")

# Per-chromosome effect size (median delta)
chr_effect <- chr_delta[sample %in% c("XG_150", "XG_151"),
                        .(sample, chrom, delta, n_intervals)]
message("Per-chromosome deltas (XG_150):")
chr150 <- chr_effect[sample == "XG_150"][order(delta)]
for (i in seq_len(min(nrow(chr150), 10))) {
  message(sprintf("  %s: delta = %.3f", chr150$chrom[i], chr150$delta[i]))
}
message("  ...")
for (i in seq(max(1, nrow(chr150) - 4), nrow(chr150))) {
  message(sprintf("  %s: delta = %.3f", chr150$chrom[i], chr150$delta[i]))
}

# Overall effect size (Cliff's delta-like: median pairwise difference)
message("\nOverall effect (median within-chromosome delta):")
for (s in c("XG_150", "XG_151")) {
  d <- chr_delta[sample == s & !is.na(delta), delta]
  message(sprintf("  %s: median = %.3f, IQR = [%.3f, %.3f]",
                  SAMPLE_LABELS[s], median(d), quantile(d, 0.25), quantile(d, 0.75)))
}

# ============================================================================
# SECTION 5: Per-interval enrichment
# ============================================================================
message("\n=== SECTION 5: Per-interval enrichment ===")

# Per-interval delta: foreground log2_signal vs background median per iteration
bg_interval_median <- bg[, .(median_bg_log2 = median(log2_signal, na.rm = TRUE)),
                           by = .(sample, interval_id)]

interval_delta <- merge(fg[, .(sample, chrom, interval_id, size, log2_signal, sample_label)],
                        bg_interval_median, by = c("sample", "interval_id"))
interval_delta[, delta := log2_signal - median_bg_log2]

# Summary
for (s in c("XG_150", "XG_151")) {
  sub <- interval_delta[sample == s]
  message(sprintf("\n%s:", SAMPLE_LABELS[s]))
  message(sprintf("  Intervals with positive delta: %d / %d (%.1f%%)",
                  sum(sub$delta > 0), nrow(sub), mean(sub$delta > 0) * 100))
  message(sprintf("  Median per-interval delta: %.3f", median(sub$delta)))
}

# ============================================================================
# SECTION 6: Figures
# ============================================================================
message("\n=== SECTION 6: Figures ===")

# --- Panel A: Chromosome-level enrichment (beeswarm + boxplot) ---
message("Panel A: Chromosome-level enrichment")

chr_delta_plot <- chr_delta[sample %in% c("XG_150", "XG_151")]
chr_delta_plot[, sample_label := factor(sample_label, levels = c("CENP-A rep1", "CENP-A rep2"))]

p_chr <- ggplot(chr_delta_plot, aes(x = "All chromosomes", y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_jitter(aes(color = delta > 0), size = 2.5, width = 0.2, alpha = 0.8) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.3, linewidth = 0.4) +
  facet_wrap(~ sample_label, nrow = 1) +
  scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"), guide = "none") +
  labs(x = NULL,
       y = expression(Delta * " median log"[2] * "(normalized signal) [349bp - shuffled bg]"),
       title = "CENP-A enrichment at 349-bp repeat loci",
       subtitle = sprintf("%d/%d chromosomes positive (Wilcoxon P = %.2e rep1, %.2e rep2)",
                          sum(chr_delta_wide$XG_150 > 0, na.rm = TRUE), n_chr,
                          wilcox_150$p.value, wilcox_151$p.value)) +
  theme_349 +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(PLOTS_DIR, "panelA_chromosome_enrichment.pdf"), p_chr, width = 8, height = 4.5)
ggsave(file.path(PLOTS_DIR, "panelA_chromosome_enrichment.png"), p_chr, width = 8, height = 4.5, dpi = 300)

# --- Panel B: Chromosome barplot ---
message("Panel B: Chromosome barplot")

p_bar <- ggplot(chr_delta_plot, aes(x = chrom_ordered, y = delta, fill = sample_label)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS, name = NULL) +
  labs(x = "Chromosome (ordered by enrichment)", y = expression(Delta * " log"[2] * "(normalized signal)"),
       title = "349-bp repeat loci: CENP-A enrichment by chromosome") +
  theme_349 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(PLOTS_DIR, "panelB_chromosome_barplot.pdf"), p_bar, width = 10, height = 4.5)
ggsave(file.path(PLOTS_DIR, "panelB_chromosome_barplot.png"), p_bar, width = 10, height = 4.5, dpi = 300)

# --- Panel C: Per-interval enrichment ---
message("Panel C: Per-interval enrichment")

interval_plot <- interval_delta[sample %in% c("XG_150", "XG_151")]
interval_plot[, size_kb := size / 1000]
interval_plot[, sample_label := factor(sample_label, levels = c("CENP-A rep1", "CENP-A rep2"))]

p_interval <- ggplot(interval_plot, aes(x = log10(size_kb), y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = chrom), size = 0.8, alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8, span = 0.8) +
  facet_wrap(~ sample_label, nrow = 1) +
  scale_color_discrete(guide = "none") +
  labs(x = expression(log[10] * "(interval size / kb)"),
       y = expression(Delta * " log"[2] * "(normalized signal) [per-interval fg - bg]"),
       title = "349-bp repeat loci: per-interval CENP-A enrichment vs interval size") +
  theme_349

ggsave(file.path(PLOTS_DIR, "panelC_per_interval_enrichment.pdf"), p_interval, width = 9, height = 4)
ggsave(file.path(PLOTS_DIR, "panelC_per_interval_enrichment.png"), p_interval, width = 9, height = 4, dpi = 300)

# --- Supp: Replicate correlation ---
message("Supplementary: Replicate correlation")

chr_wide <- dcast(chr_delta_wide, chrom ~ ., value.var = c("XG_150", "XG_151"))
names(chr_wide) <- c("chrom", "rep1", "rep2")
rho <- cor(chr_wide$rep1, chr_wide$rep2, method = "spearman", use = "complete.obs")

p_rep <- ggplot(chr_wide, aes(x = rep1, y = rep2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "#2166AC", fill = "#2166AC", alpha = 0.15, linewidth = 0.6) +
  geom_point(size = 2.5, alpha = 0.8, color = "#2166AC") +
  annotate("text", x = min(chr_wide$rep1) + 0.1, y = max(chr_wide$rep2) - 0.1,
           label = paste("Spearman rho =", round(rho, 3)),
           hjust = 0, size = 3.5) +
  labs(x = expression(Delta * " log"[2] * "(normalized signal) CENP-A rep1"),
       y = expression(Delta * " log"[2] * "(normalized signal) CENP-A rep2"),
       title = "349-bp repeat loci: CENP-A replicate concordance") +
  theme_349

ggsave(file.path(PLOTS_DIR, "supp_replicate_correlation.pdf"), p_rep, width = 5, height = 5)
ggsave(file.path(PLOTS_DIR, "supp_replicate_correlation.png"), p_rep, width = 5, height = 5, dpi = 300)

# --- Supp: H3K27ac comparison ---
message("Supplementary: H3K27ac comparison")

chr_h3k <- chr_delta[sample %in% c("XG_152", "XG_153")]
chr_h3k_wide <- dcast(chr_h3k, chrom ~ sample, value.var = "delta")
if (nrow(chr_h3k_wide) > 0 && all(c("XG_152", "XG_153") %in% names(chr_h3k_wide))) {
  p_h3k <- ggplot(chr_h3k_wide, aes(x = XG_152, y = XG_153)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", size = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", size = 0.3) +
    geom_point(size = 2.5, alpha = 0.8, color = "#B2182B") +
    labs(x = expression(Delta * " log"[2] * "(normalized signal) H3K27ac"),
         y = expression(Delta * " log"[2] * "(normalized signal) H3K27ac rep2"),
         title = "349-bp repeat loci: H3K27ac per chromosome") +
    theme_349

  ggsave(file.path(PLOTS_DIR, "supp_h3k27ac_comparison.pdf"), p_h3k, width = 5, height = 5)
  ggsave(file.path(PLOTS_DIR, "supp_h3k27ac_comparison.png"), p_h3k, width = 5, height = 5, dpi = 300)
}

# --- Supp: Interval size distribution ---
message("Supplementary: Interval size distribution")

p_size <- ggplot(fg[sample == "XG_150"], aes(x = size / 1000)) +
  geom_histogram(bins = 60, fill = "#2166AC", alpha = 0.7, color = "white", linewidth = 0.1) +
  scale_x_log10() +
  labs(x = "Interval size (kb, log scale)", y = "Count",
       title = "349-bp repeat loci: merged interval size distribution") +
  theme_349

ggsave(file.path(PLOTS_DIR, "supp_interval_sizes.pdf"), p_size, width = 6, height = 4)
ggsave(file.path(PLOTS_DIR, "supp_interval_sizes.png"), p_size, width = 6, height = 4, dpi = 300)

# --- Supp: Per-chromosome interval counts ---
message("Supplementary: Per-chromosome interval counts")

chr_counts <- fg[sample == "XG_150", .N, by = chrom]
chr_counts[, chrom_factor := factor(chrom, levels = rev(CHROMOSOMES))]
p_counts <- ggplot(chr_counts, aes(x = chrom_factor, y = N)) +
  geom_col(fill = "#2166AC", alpha = 0.8) +
  coord_flip() +
  labs(x = NULL, y = "Number of 349-bp intervals",
       title = "349-bp repeat loci per chromosome") +
  theme_349

ggsave(file.path(PLOTS_DIR, "supp_interval_counts.pdf"), p_counts, width = 5, height = 6)
ggsave(file.path(PLOTS_DIR, "supp_interval_counts.png"), p_counts, width = 5, height = 6, dpi = 300)

# ============================================================================
# SECTION 7: Save results
# ============================================================================
message("\n=== SECTION 7: Saving results ===")

# Save summary
sink(file.path(RESULTS_DIR, "349bp_results_summary.txt"))
cat("CENP-A CUT&Tag Enrichment at 349-bp Repeat Loci — Results Summary\n")
cat("=================================================================\n\n")
cat(sprintf("Date: %s\n\n", Sys.time()))
cat(sprintf("Samples: %d (XG_150 CENP-A rep1, XG_151 CENP-A rep2, XG_152 H3K27ac, XG_153 H3K27ac rep2)\n", length(SAMPLES)))
cat(sprintf("Pseudocount: %.2e\n\n", pseudocount))
cat(sprintf("Merged 349-bp intervals: %d\n", uniqueN(fg$interval_id)))
cat(sprintf("Background: %d iterations, %d total shuffled intervals\n\n", n_iter, nrow(bg)))

cat("Primary inference (chromosome-level, fg vs shuffled bg):\n")
cat(sprintf("  XG_150: V = %.0f, P = %.3e\n", wilcox_150$statistic, wilcox_150$p.value))
cat(sprintf("  XG_151: V = %.0f, P = %.3e\n", wilcox_151$statistic, wilcox_151$p.value))
cat(sprintf("  Positive chromosomes (XG_150): %d / %d (%.1f%%)\n",
            sum(chr_delta_wide$XG_150 > 0, na.rm = TRUE), n_chr, prop_pos_150 * 100))
cat(sprintf("  Positive chromosomes (XG_151): %d / %d (%.1f%%)\n",
            sum(chr_delta_wide$XG_151 > 0, na.rm = TRUE), n_chr, prop_pos_151 * 100))

cat("\nPer-interval enrichment (fg vs bg median):\n")
for (s in c("XG_150", "XG_151")) {
  sub <- interval_delta[sample == s]
  cat(sprintf("  %s: %d / %d intervals positive (%.1f%%), median delta = %.3f\n",
              SAMPLE_LABELS[s], sum(sub$delta > 0, na.rm = TRUE), nrow(sub),
              mean(sub$delta > 0, na.rm = TRUE) * 100, median(sub$delta, na.rm = TRUE)))
}

cat(sprintf("\nCENP-A replicate Spearman rho: %.3f\n", rho))
sink()

message("Summary saved to: ", file.path(RESULTS_DIR, "349bp_results_summary.txt"))

# Per-chromosome statistics
stats_out <- chr_delta_wide
stats_out$n_intervals <- chr_fg$n_intervals[match(stats_out$chrom, chr_fg$chrom[chr_fg$sample == "XG_150"])]
names(stats_out) <- c("chromosome", "XG_150_delta", "XG_151_delta", "n_intervals")
stats_out <- stats_out[order(stats_out$XG_150_delta), ]
fwrite(stats_out, file.path(RESULTS_DIR, "349bp_statistics.csv"))
message("Statistics saved to: ", file.path(RESULTS_DIR, "349bp_statistics.csv"))

# Session info
writeLines(capture.output(sessionInfo()), file.path(RESULTS_DIR, "session_info.txt"))

message("\n=== Analysis complete ===")
message("Plots: ", PLOTS_DIR)
message("Results: ", RESULTS_DIR)
