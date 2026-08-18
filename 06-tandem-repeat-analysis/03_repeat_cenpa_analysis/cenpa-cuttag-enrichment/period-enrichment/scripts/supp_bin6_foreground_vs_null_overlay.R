#!/usr/bin/env Rscript
# ============================================================================
# supp_bin6_foreground_vs_null_overlay.R
# PLOTTING half of the per-interval overlay for the 348-349 bp monomer
# (bin 6):
#   * grey  = NULL per-interval signal  (all shuffled intervals across 1,000
#             iterations, pooled)
#   * color = observed CENP-A foreground per-interval signal
#   Both shown as densities (area = 1), so the different N (900 foreground vs
#   ~840k null) does not distort the comparison. Median vertical lines mark each
#   distribution's location; P_adj from the 9-bin permutation results is annotated.
#
# This script ONLY renders the figure. All the expensive computation (reading
# the per-interval signal + bin-6 permutation null files, pooling the 1,000
# shuffled iterations) is done by
# scripts/supp_bin6_foreground_vs_null_overlay_prepare.R, which caches the
# per-interval log2 signals in
# ../results/supp_bin6_foreground_vs_null_distribution_data.tsv.
#
# Reading that table and drawing the densities takes seconds, so you can
# iterate freely on plot format (colors, labels, sizes, layout) without
# re-reading the raw input files.
#
# Normalization matches 3_analyze_period_enrichment.R EXACTLY (see the
# _prepare.R header for the formula).
#
# Usage:
#   conda activate r-visualizations
#   Rscript supp_bin6_foreground_vs_null_overlay.R
#
#   (If the cache is missing, run the _prepare.R script first.)
#
# Saved: ../plots/supp_bin6_foreground_vs_null_distribution.{pdf,png}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

PERIOD_DIR  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PLOTS_DIR   <- file.path(PERIOD_DIR, "plots")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

SAMPLES <- c("XG_150", "XG_151")                       # CENP-A replicates
SAMPLE_LABELS <- c("XG_150" = "CENP-A rep1", "XG_151" = "CENP-A rep2")
BIN6 <- 6  # 348-349 bp = the 349-bp monomer family in the 9-bin scheme

DATA_FILE <- file.path(RESULTS_DIR,
                       "supp_bin6_foreground_vs_null_distribution_data.tsv")
PERM_CSV  <- file.path(RESULTS_DIR,
                       "period_enrichment_permutation_results.csv")

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

# ── Read the cached per-interval values (written by ..._prepare.R) ────────────
if (!file.exists(DATA_FILE)) {
  stop("Cached plotting data not found: ", DATA_FILE,
       "\n  Run scripts/supp_bin6_foreground_vs_null_overlay_prepare.R first ",
       "(it reads the signal/permutation files and writes this table).")
}
d <- fread(DATA_FILE, sep = "\t")
d[, sample_label := factor(SAMPLE_LABELS[sample], levels = c("CENP-A rep1", "CENP-A rep2"))]
fg_bin6 <- d[grp == "foreground"]
nul_bin6 <- d[grp == "null"]

message(sprintf("Foreground bin-6 intervals/sample: %d", nrow(fg_bin6) / length(SAMPLES)))
message(sprintf("Null bin-6 rows (interval x iteration) per sample: %d",
                nrow(nul_bin6) / length(SAMPLES)))

# ── Medians for the vertical markers ───────────────────────────────────────
fg_med <- fg_bin6[, .(median = median(log2_signal, na.rm = TRUE)), by = .(sample, sample_label)]
nul_med <- nul_bin6[, .(median = median(log2_signal, na.rm = TRUE)), by = .(sample, sample_label)]
message("Foreground medians:")
print(fg_med)
message("Null medians:")
print(nul_med)

# ── P-values from the 9-bin permutation results (small pre-computed CSV) ────
perm_pvals <- fread(PERM_CSV)
pvals_bin6 <- perm_pvals[bin_id == BIN6 & sample %in% SAMPLES]
pvals_bin6[, sample_label := factor(SAMPLE_LABELS[sample], levels = c("CENP-A rep1", "CENP-A rep2"))]

# ── Figure: per-interval density overlay ───────────────────────────────────
p_overlay <- ggplot() +
  # Null per-interval distribution (reference)
  geom_density(data = nul_bin6,
               aes(x = log2_signal, fill = "Null (shuffled)"),
               color = NA, alpha = 0.7, linewidth = 0.3) +
  # Observed foreground per-interval distribution
  geom_density(data = fg_bin6,
               aes(x = log2_signal, fill = sample_label),
               color = "grey20", linewidth = 0.25, alpha = 0.55) +
  # Median markers
  geom_vline(data = nul_med, aes(xintercept = median),
             color = "grey30", linetype = "dashed", linewidth = 0.4) +
  geom_vline(data = fg_med, aes(xintercept = median, color = sample_label),
             linetype = "dashed", linewidth = 0.6) +
  # P-value annotation (depletion: observed < null)
  geom_text(data = pvals_bin6,
            aes(x = -2.6, y = Inf,
                label = sprintf("P_adj = %.3f %s", p_adj, sig)),
            hjust = 0, vjust = 1.3, size = 3, fontface = "bold") +
  facet_wrap(~ sample_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(name = NULL,
                    values = c("Null (shuffled)" = "grey60",
                               "CENP-A rep1" = "#2166AC",
                               "CENP-A rep2" = "#92C5DE")) +
  scale_color_manual(values = c("CENP-A rep1" = "#2166AC",
                                "CENP-A rep2" = "#92C5DE"),
                     guide = "none") +
  coord_cartesian(xlim = c(-10.5, 0.5)) +
  labs(x = expression("Per-interval log"[2] * "(normalized signal)"),
       y = "Density",
       title = "348-349 bp monomer (bin 6): observed vs null per-interval distribution",
       subtitle = paste0("Grey: null (all shuffled intervals, 1,000 iterations pooled). ",
                         "Colored: observed CENP-A intervals (n = 900).\n",
                         "Dashed: medians. Spike at log2 = -9.4 = zero-coverage intervals. ",
                         "Observed median < null median: depletion.")) +
  theme_period

ggsave(file.path(PLOTS_DIR, "supp_bin6_foreground_vs_null_distribution.pdf"),
       p_overlay, width = 7, height = 6)
ggsave(file.path(PLOTS_DIR, "supp_bin6_foreground_vs_null_distribution.png"),
       p_overlay, width = 7, height = 6, dpi = 300)
message("Saved: supp_bin6_foreground_vs_null_distribution.pdf/png")
message("Done.")
