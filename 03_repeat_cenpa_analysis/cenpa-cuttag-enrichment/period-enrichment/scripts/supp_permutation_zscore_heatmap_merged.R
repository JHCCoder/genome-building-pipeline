#!/usr/bin/env Rscript
# ============================================================================
# supp_permutation_zscore_heatmap_merged.R
# PLOTTING half of the per-chromosome x per-bin CENP-A enrichment Z-score
# heatmap (MERGED arrays, rep1 = XG_150).
#
# This script ONLY renders the figure. All the expensive computation
# (per-(chrom,bin) Z-scores + per-chromosome permutation significance) is done
# by scripts/supp_permutation_zscore_heatmap_merged_prepare.R, which caches the
# result in ../results/supp_permutation_zscore_heatmap_merged_matrix.tsv.
#
# Reading that small table and drawing the heatmap takes seconds, so you can
# iterate freely on plot format (colors, labels, sizes, ordering) without
# re-running the ~GB permutation input files.
#
# Paper version changes (kept in the plot):
#
#   1. Microsatellite (bin 1) and minisatellite (bin 2) COLUMNS REMOVED —
#      only bins 3-9 are plotted.
#   2. Significance stars INSIDE THE CELL for Z-scores that are significant by
#      a PER-CHROMOSOME permutation test (same empirical test as the bin-level
#      analysis but restricted to each chromosome; matched design for ctrl bins
#      3-5, full-set for bins 6-9; BH-FDR within each bin, rep1).
#   3. Rows sorted CHR1 -> CHR28 -> CHRX (natural order), CHR1 AT THE TOP,
#      not by bin-7 z-score.
#   4. Missing cells (e.g. chromosomes with zero 348-349 bp arrays) are drawn
#      as explicit grey NA tiles instead of being silently absent.
#   5. X-axis (period bin) labels at the BOTTOM (not the top).
#   6. Color scale flipped: RED = enriched (Z>0), BLUE = depleted (Z<0).
#   7. Short title, NO subtitle — the explanatory text lives in a notebook
#      markdown cell (figure_notebook.ipynb, heatmap caption).
#
# For the methodology (Z-score definition, per-chromosome test, why some
# bin-6 cells are grey) see the header of the _prepare.R script.
#
# Usage:
#   conda activate r-visualizations
#   Rscript supp_permutation_zscore_heatmap_merged.R             # default 11 x 8 in
#   Rscript supp_permutation_zscore_heatmap_merged.R 12 9        # custom size (inches)
#
#   (If the cached matrix is missing, run the _prepare.R script first.)
#
# Saved: ../plots/supp_permutation_zscore_heatmap_merged.{pdf,png,svg}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

PERIOD_DIR  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PLOTS_DIR   <- file.path(PERIOD_DIR, "plots")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")
OUT_SUFFIX  <- "merged"

MATRIX_FILE <- file.path(RESULTS_DIR,
                         "supp_permutation_zscore_heatmap_merged_matrix.tsv")

# ── Figure dimensions (inches) — single knobs; edit + re-run to resize ────────
# Optional command-line override (inches):
#   Rscript supp_permutation_zscore_heatmap_merged.R WIDTH HEIGHT
DEFAULT_WIDTH  <- 11
DEFAULT_HEIGHT <- 8

args <- commandArgs(trailingOnly = TRUE)
PLOT_WIDTH  <- if (length(args) >= 1) as.numeric(args[1]) else DEFAULT_WIDTH
PLOT_HEIGHT <- if (length(args) >= 2) as.numeric(args[2]) else DEFAULT_HEIGHT
if (is.na(PLOT_WIDTH) || is.na(PLOT_HEIGHT) ||
    PLOT_WIDTH <= 0 || PLOT_HEIGHT <= 0) {
  stop("Invalid WIDTH/HEIGHT. Usage: Rscript supp_permutation_zscore_heatmap_merged.R WIDTH HEIGHT",
       "\n  WIDTH and HEIGHT are positive numbers in inches.")
}

BIN_SHORT <- c(
  "3" = "51-192 bp", "4" = "193-195 bp", "5" = "196-347 bp",
  "6" = "348-349 bp", "7" = "350-385 bp", "8" = "386-390 bp", "9" = "391+ bp"
)

# ── Read the cached Z-score matrix (written by ..._prepare.R) ─────────────────
if (!file.exists(MATRIX_FILE)) {
  stop("Cached Z-score matrix not found: ", MATRIX_FILE,
       "\n  Run scripts/supp_permutation_zscore_heatmap_merged_prepare.R first ",
       "(it reads the raw signal/permutation files and writes this table).")
}
chr_z <- fread(MATRIX_FILE, sep = "\t")
PLOT_BINS <- sort(unique(chr_z$bin_id))

# ── Rebuild view factors from the cached data ─────────────────────────────────
# Natural chromosome order: chr1..chr28, chrX (chrY if it had data). chr1 at TOP.
chr_set <- unique(chr_z$chrom)
num <- suppressWarnings(as.numeric(sub("chr", "", chr_set)))
num[grepl("X", chr_set)] <- 100
num[grepl("Y", chr_set)] <- 101
chr_order <- chr_set[order(num)]
chr_z[, chrom_factor := factor(chrom, levels = rev(chr_order))]
chr_z[, bin_factor := factor(bin_id, levels = PLOT_BINS)]

# Bin labels (no stars on the axis; stars are inside cells)
x_labels <- setNames(BIN_SHORT[as.character(PLOT_BINS)], as.character(PLOT_BINS))

# ── Plot ──────────────────────────────────────────────────────────────────────
p_zheat <- ggplot(chr_z, aes(x = bin_factor, y = chrom_factor, fill = z_capped)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(z_score), "",
                               paste0(sprintf("%.1f", z_score), sig_mark))),
            size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, na.value = "grey90",
                       name = "CENP-A signal\nZ-score vs\npermutation null") +
  scale_x_discrete(labels = x_labels, position = "bottom",
                   expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "TRF period bin",
       y = "Chromosome",
       title = "Chromosome-specific CENP-A signal by repeat-period bin") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 10),
        strip.background = element_rect(fill = "grey95", color = "grey80"),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black"),
            plot.title = element_text(size = 12, face = "bold", margin = margin(b = 4)),
	axis.title = element_text(color = "black"))

ggsave(file.path(PLOTS_DIR, paste0("supp_permutation_zscore_heatmap_", OUT_SUFFIX, ".pdf")),
       p_zheat, width = PLOT_WIDTH, height = PLOT_HEIGHT)
ggsave(file.path(PLOTS_DIR, paste0("supp_permutation_zscore_heatmap_", OUT_SUFFIX, ".png")),
       p_zheat, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)
ggsave(file.path(PLOTS_DIR, paste0("supp_permutation_zscore_heatmap_", OUT_SUFFIX, ".svg")),
       p_zheat, width = PLOT_WIDTH, height = PLOT_HEIGHT, device = svglite::svglite)
message(sprintf("Saved: supp_permutation_zscore_heatmap_%s.pdf/png/svg (%g x %g in)",
                OUT_SUFFIX, PLOT_WIDTH, PLOT_HEIGHT))
message("Done.")
