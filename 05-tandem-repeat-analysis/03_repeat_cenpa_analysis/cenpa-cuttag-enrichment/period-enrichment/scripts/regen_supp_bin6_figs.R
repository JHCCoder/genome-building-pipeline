#!/usr/bin/env Rscript
# ============================================================================
# regen_supp_bin6_figs.R — ONE-OFF correction of the period-enrichment
# supplementary figures after the 2026-08-02 rebin to 9 bins.
#
# In the 9-bin scheme the 349-bp monomer family is bin 6 (348-349 bp), NOT
# bin 7 (350-385 bp). The figures named "349-bp monomer (bin 7)" were plotting
# the wrong bin. 3_analyze_period_enrichment.R has been fixed to use bin 6 for
# future full runs; this script regenerates the affected figures NOW without
# re-running the whole (heavy) analysis:
#   * supp_bin6_null_distribution.{pdf,png}  (was supp_bin7_null_distribution)
#   * supp_349bp_copies_vs_signal.{pdf,png}  (now bin 6)
#   * supp_349bp_match_vs_signal.{pdf,png}   (now bin 6)
#
# Usage:
#   conda activate r-visualizations
#   Rscript regen_supp_bin6_figs.R
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

PERIOD_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PARENT_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
COUNTS_DIR <- file.path(PERIOD_DIR, "data", "counts")
PER_BIN_DIR <- file.path(PERIOD_DIR, "data", "permutation", "per_bin")
PLOTS_DIR  <- file.path(PERIOD_DIR, "plots")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

SAMPLES <- c("XG_150", "XG_151", "XG_152", "XG_153")
CENPA <- c("XG_150", "XG_151")
SAMPLE_LABELS <- c("XG_150"="CENP-A rep1", "XG_151"="CENP-A rep2",
                   "XG_152"="H3K27ac", "XG_153"="H3K27ac rep2")
SAMPLE_COLORS <- c("XG_150"="#2166AC", "XG_151"="#92C5DE",
                   "XG_152"="#B2182B", "XG_153"="#D6604D")
BIN6 <- 6  # 348-349 bp = the 349-bp monomer family in the 9-bin scheme

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

# ── Library sizes ───────────────────────────────────────────────────────────
lib_file <- file.path(PARENT_DIR, "data", "qc", "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
} else {
  lib_sizes <- c("XG_150"=27156788, "XG_151"=25305681, "XG_152"=25305681, "XG_153"=25305681)
}

# ── Foreground: load signal for ALL samples (pseudo must match the analysis) ─
signal_list <- list()
for (s in SAMPLES) {
  x <- fread(file.path(COUNTS_DIR, paste0("trf_signal_", s, ".tsv")), header = FALSE,
             col.names = c("chrom","start","end","interval_id","bin_id","period_size",
                           "copies_aligned","match_percent","mean_coverage"))
  x[, sample := s]
  signal_list[[s]] <- x
}
signal <- rbindlist(signal_list)
signal[, lib_size := lib_sizes[sample]]
signal[, norm_signal := mean_coverage / (lib_size / 1e6)]
pseudo <- min(signal$norm_signal[signal$norm_signal > 0], na.rm = TRUE) / 2
signal[, log2_signal := log2(norm_signal + pseudo)]
message(sprintf("Pseudocount: %.4e", pseudo))

# Bin-6 foreground median per CENP-A sample
fg_bin6 <- signal[bin_id == BIN6 & sample %in% CENPA,
                  .(fg_median = median(log2_signal, na.rm = TRUE)), by = .(sample, bin_id)]
fg_bin6[, sample_label := factor(SAMPLE_LABELS[sample],
                                 levels = c("CENP-A rep1", "CENP-A rep2"))]

# Bin-6 per-iteration null medians (full-set null, from per-bin files)
bg_bin6 <- rbindlist(lapply(CENPA, function(s) {
  f <- file.path(PER_BIN_DIR, "full", sprintf("bin%d_%s.tsv", BIN6, s))
  x <- fread(f, header = FALSE,
             col.names = c("chrom","start","end","bin_id","iter","mean_coverage"))
  x[, sample := s]
  x[, norm_signal := mean_coverage / (lib_sizes[s] / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
  x[, .(bg_median = median(log2_signal, na.rm = TRUE)), by = .(sample, bin_id, iter)]
}))
bg_bin6[, sample_label := factor(SAMPLE_LABELS[sample],
                                 levels = c("CENP-A rep1", "CENP-A rep2"))]

# P-values from the (correct) 9-bin results CSV
perm_pvals <- fread(file.path(RESULTS_DIR, "period_enrichment_permutation_results.csv"))

# ══ Supp B: bin-6 null distribution vs observed CENP-A ═════════════════════
p_null <- ggplot(bg_bin6[sample %in% CENPA]) +
  geom_histogram(aes(x = bg_median, fill = sample_label),
                 bins = 60, alpha = 0.6, position = "identity") +
  geom_vline(data = fg_bin6[sample %in% CENPA],
             aes(xintercept = fg_median, color = sample_label),
             linewidth = 1.2) +
  geom_text(data = perm_pvals[bin_id == BIN6 & sample %in% CENPA],
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

ggsave(file.path(PLOTS_DIR, "supp_bin6_null_distribution.pdf"), p_null, width = 7, height = 6)
ggsave(file.path(PLOTS_DIR, "supp_bin6_null_distribution.png"), p_null, width = 7, height = 6, dpi = 300)
message("Saved: supp_bin6_null_distribution.pdf/png")

# ══ Supp C / D: copy number & repeat homogeneity vs CENP-A (bin 6) ═════════
bin6_signal <- signal[bin_id == BIN6 & sample %in% CENPA]
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

message("Done. Bin-6 supp figures regenerated.")
