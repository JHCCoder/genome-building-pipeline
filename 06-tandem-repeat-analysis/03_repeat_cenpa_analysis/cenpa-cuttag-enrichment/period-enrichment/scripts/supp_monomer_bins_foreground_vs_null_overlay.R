#!/usr/bin/env Rscript
# ============================================================================
# supp_monomer_bins_foreground_vs_null_overlay.R
# PLOTTING half of the per-interval overlay for the three CENP-A monomer
# families: bins 4 (193-195 bp), 6 (348-349 bp) and 8 (386-390 bp).
#
# This script ONLY renders the figure. All the expensive computation (reading
# the per-interval signal + permutation null files, pooling the 1,000 shuffled
# iterations, thinning nulls >1M rows) is done by
# scripts/supp_monomer_bins_foreground_vs_null_overlay_prepare.R, which caches
# the per-interval values in
# ../results/supp_monomer_bins_foreground_vs_null_distribution_data.tsv.
#
# Reading that table and drawing the densities takes seconds, so you can
# iterate freely on plot format (colors, labels, sizes, layout) without
# re-reading the raw input files.
#
# Figure content (unchanged):
#   * grey  = NULL per-interval signal (all shuffled intervals across 1,000
#             iterations, pooled)
#   * color = observed CENP-A foreground per-interval signal
#   Both shown as densities (area = 1) in a 3 (bin) x 2 (replicate) panel grid.
#   A dashed colored line marks the OBSERVED median per panel (descriptive only;
#   the grey null median is deliberately NOT drawn - it would invite the false
#   reading that the observed median was tested against that one line). The panel
#   label (enriched / depleted) is taken from the PERMUTATION test (sign of the
#   matched effect), so it always agrees with the reported P-value.
#
# Normalization matches 3_analyze_period_enrichment.R EXACTLY (see the
# _prepare.R header for the formula). NOTE on the zero-coverage spike at
# log2 = -9.38: a real property of the shuffled background, not an artifact.
#
# Usage:
#   conda activate r-visualizations
#   Rscript supp_monomer_bins_foreground_vs_null_overlay.R
#
#   (If the cache is missing, run the _prepare.R script first.)
#
# Saved: ../plots/supp_monomer_bins_foreground_vs_null_distribution.{pdf,png,svg}
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
BINS <- c(4, 6, 8)   # 193-195 bp, 348-349 bp, 386-390 bp
BIN_LABELS <- c("4" = "193-195 bp (bin 4)", "6" = "348-349 bp (bin 6)", "8" = "386-390 bp (bin 8)")

DATA_FILE <- file.path(RESULTS_DIR,
                       "supp_monomer_bins_foreground_vs_null_distribution_data.tsv")
PERM_CSV  <- file.path(RESULTS_DIR,
                       "period_enrichment_permutation_results.csv")

# ── Read the cached per-interval values (written by ..._prepare.R) ────────────
if (!file.exists(DATA_FILE)) {
  stop("Cached plotting data not found: ", DATA_FILE,
       "\n  Run scripts/supp_monomer_bins_foreground_vs_null_overlay_prepare.R first ",
       "(it reads the signal/permutation files and writes this table).")
}
d <- fread(DATA_FILE, sep = "\t")
d[, bin_factor := factor(BIN_LABELS[as.character(bin_id)], levels = BIN_LABELS[as.character(BINS)])]
d[, sample_label := factor(SAMPLE_LABELS[sample], levels = c("CENP-A rep1", "CENP-A rep2"))]

# ── Medians (descriptive; only the OBSERVED median line is drawn) ──────────
med <- d[, .(med = median(value, na.rm = TRUE)), by = .(bin_id, bin_factor, sample, sample_label, grp)]
med_wide <- dcast(med, bin_id + bin_factor + sample + sample_label ~ grp, value.var = "med")

# ── Direction from the PERMUTATION test (sign of the matched effect D) ─────
# Small pre-computed results CSV (not raw data) — read here, not in _prepare,
# so the plot script stays self-contained. The on-panel label can never
# contradict the reported P-value.
perm_pvals <- fread(PERM_CSV)
dir <- merge(med_wide,
             perm_pvals[bin_id %in% BINS & sample %in% SAMPLES,
                        .(bin_id, sample, direction = ifelse(paired_effect > 0, "↑ Enriched", "↓ Depleted"))],
             by = c("bin_id", "sample"))
n_sig <- perm_pvals[bin_id %in% BINS & sample %in% SAMPLES & p_adj < 0.001, .N]

message("Observed vs null medians (log2 normalized; descriptive) and direction (from permutation test):")
print(copy(dir)[, .(bin_id, sample_label, direction, foreground = round(foreground, 3), null = round(null, 3))])
message(sprintf("Panels with P_adj < 0.001: %d / %d", n_sig, length(BINS) * length(SAMPLES)))

# ── Combined figure: 3 (bin) x 2 (sample) panel grid ───────────────────────
p_overlay <- ggplot(d) +
  geom_density(data = d[grp == "null"], aes(x = value, fill = "Null (shuffled)"),
               color = NA, alpha = 0.75) +
  geom_density(data = d[grp == "foreground"], aes(x = value, fill = sample_label),
               color = "grey20", linewidth = 0.25, alpha = 0.55) +
  # Observed median only. The grey null median is intentionally NOT drawn:
  # it is a pooled descriptive summary, not the statistic the permutation
  # test compares the observed median against (that is the per-iteration
  # matched difference D_j = median(observed_j) - median(shuffled_j)).
  geom_vline(data = med[grp == "foreground"], aes(xintercept = med, color = sample_label),
             linetype = "dashed", linewidth = 0.55) +
  geom_text(data = dir, aes(x = -2.4, y = Inf, label = direction),
            hjust = 0, vjust = 1.4, size = 3, fontface = "bold", color = "grey20") +
  facet_grid(bin_factor ~ sample_label, scales = "free_y") +
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
       title = "CENP-A monomer families: observed vs null per-interval distribution",
       subtitle = paste0("Grey: interval-level signal at chromosome- and length-matched shuffled placements\n",
                         "(1,000 iterations, pooled). Colored: observed intervals; dashed colored line =\n",
                         "observed median (descriptive). Spike at log2 = -9.4 = zero-coverage intervals.\n",
                         "Direction labels and P_adj < 0.001 are from the permutation test of iteration-specific\n",
                         "matched observed-vs-shuffled median differences - not from these curves.")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
        strip.background = element_rect(fill = "grey95", color = "grey80"),
        strip.text = element_text(size = 9, face = "bold"),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave(file.path(PLOTS_DIR, "supp_monomer_bins_foreground_vs_null_distribution.pdf"),
       p_overlay, width = 9, height = 7)
ggsave(file.path(PLOTS_DIR, "supp_monomer_bins_foreground_vs_null_distribution.png"),
       p_overlay, width = 9, height = 7, dpi = 300)
ggsave(file.path(PLOTS_DIR, "supp_monomer_bins_foreground_vs_null_distribution.svg"),
       p_overlay, width = 9, height = 7, device = svglite::svglite)
message("Saved: supp_monomer_bins_foreground_vs_null_distribution.pdf/png/svg")
message("Done.")
