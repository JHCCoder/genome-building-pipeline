#!/usr/bin/env Rscript
# ============================================================================
# period_violin_array.R — PLOTTING half of the horizontal per-merged-ARRAY
#                    Δ signal violin: observed log2 signal centered on the
#                    permutation null baseline, by TRF period bin.
#
#   DESCRIPTIVE layer. Each array's centered value is
#     Δ_i = log2(CPM_i + pseudo) − null_median(bin, sample)
#   where null_median = mean of the 1,000 per-iteration shuffled medians
#   (chr- and length-matched), computed on MERGED REPEAT ARRAYS (bedtools
#   merge -d 0 per bin) — the same unit as the observed values:
#     Bins 1–5: 1,000 iterations on a 5,000-array subsample
#     Bins 6–9: 1,000 iterations on all arrays (full set)
#   Each point = one merged repeat array. White dot = median.
#   Significance stars are the FORMAL layer — empirical permutation p-values
#   read from the merged-array results CSV.
#
# This script ONLY renders the figure. All the expensive computation (reading
# the merged-array signal + ~1 GB permutation background files, computing the
# per-bin null medians and per-array Δ, and the significance-star positions) is
# done by scripts/period_violin_array_prepare.R, which caches:
#     ../results/period_violin_<group>_merged_delta.tsv   (array_id, bin_id, sample, delta)
#     ../results/period_violin_<group>_merged_stars.tsv   (bin_id, sample, label, x, y)
#
# Reading those tables and drawing the violin takes seconds, so you can iterate
# freely on plot format (colors, labels, sizes, theme). NOTE: star x/y
# positions are cached from the violin geometry — if you change DODGE_W or any
# layout that shifts the violins, re-run the _prepare script to re-align them.
# FONT_SIZE changes never need a re-run.
#
# Usage:
#   conda activate r-visualizations
#   Rscript period_violin_array.R CENP-A    # writes period_violin_cenpa_merged.{pdf,png,svg}
#   Rscript period_violin_array.R H3K27ac   # writes period_violin_h3k27ac_merged.{pdf,png,svg}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

PERIOD_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

# ── Sample group (CENP-A default; H3K27ac for the H3K27ac violin) ────────────
args <- commandArgs(trailingOnly = TRUE)
GROUP <- if (length(args) >= 1) args[1] else "CENP-A"
if (!GROUP %in% c("CENP-A", "H3K27ac")) stop("Unknown GROUP: ", GROUP)

GROUP_CONFIG <- list(
  "CENP-A"  = list(samples = c("XG_150", "XG_151"), label = "CENP-A",
                   colors = c("#1B4F72", "#7FB3D8"), prefix = "period_violin_cenpa_merged"),
  "H3K27ac" = list(samples = c("XG_152", "XG_153"), label = "H3K27ac",
                   colors = c("#B2182B", "#F4A582"), prefix = "period_violin_h3k27ac_merged")
)
cfg         <- GROUP_CONFIG[[GROUP]]
SAMPLES     <- cfg$samples
GROUP_LABEL <- cfg$label
GROUP_COLORS <- setNames(cfg$colors, paste0(GROUP_LABEL, "  rep", 1:2))
OUT_PREFIX  <- cfg$prefix

CENTROMERIC_BINS  <- 6:9
CONTROL_BINS      <- 1:5
ALL_BINS          <- c(CONTROL_BINS, CENTROMERIC_BINS)
DODGE_W           <- 0.8

# ── Font size (pt) — single knob; change this one value to scale all text ────
FONT_SIZE <- 12

BIN_LABELS <- c(
  "1" = "1–10 bp\nmicrosatellites",
  "2" = "11–50 bp\nminisatellites",
  "3" = "51–192 bp",
  "4" = "193–195 bp",
  "5" = "196–347 bp",
  "6" = "348–349 bp",
  "7" = "350–385 bp",
  "8" = "386–390 bp",
  "9" = "391+ bp"
)

DELTA_FILE <- file.path(RESULTS_DIR, paste0(cfg$prefix, "_delta.tsv"))
STARS_FILE <- file.path(RESULTS_DIR, paste0(cfg$prefix, "_stars.tsv"))

# ── Read the cached per-array Δ signal and star positions ─────────────────────
if (!file.exists(DELTA_FILE)) {
  stop("Cached Δ-signal table not found: ", DELTA_FILE,
       "\n  Run scripts/period_violin_array_prepare.R ", GROUP, " first ",
       "(it reads the signal/permutation files and writes this table).")
}
dt <- fread(DELTA_FILE, sep = "\t")
sig_stars <- fread(STARS_FILE, sep = "\t")

# ── Rebuild view factors from the cached data (same as the analysis) ──────────
# NOTE: per-bin merged-array counts are no longer shown on the y axis — they
# moved to the companion barplot strip (period_array_count_length.R).
BIN_LABELS_N <- setNames(
  vapply(ALL_BINS, function(b) {
    BIN_LABELS[as.character(b)]
  }, character(1)),
  as.character(ALL_BINS)
)

dt[, bin_label := factor(BIN_LABELS_N[as.character(bin_id)],
                          levels = BIN_LABELS_N[as.character(ALL_BINS)])]
dt[, sample_label := paste0(GROUP_LABEL, "  rep", match(sample, SAMPLES))]

bin_summary <- dt[, .(med_delta = median(delta, na.rm = TRUE), n = .N),
                  by = .(sample, bin_id, bin_label)]

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(dt, aes(x = delta, y = bin_label)) +

  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.45, color = "grey45") +

  geom_violin(
    aes(fill = sample_label),
    color = NA,
    alpha = 0.55,
    scale = "width",
    trim = TRUE,
    position = position_dodge(width = DODGE_W)
  ) +

  stat_summary(
    aes(group = sample_label),
    fun = median,
    geom = "point",
    shape = 21,
    size = 1.8,
    fill = "white",
    color = "grey20",
    stroke = 0.4,
    position = position_dodge(width = DODGE_W)
  ) +

  geom_text(data = sig_stars, aes(x = x, y = y, label = label),
            size = 0.42 * FONT_SIZE, fontface = "bold", color = "black",
            hjust = 0.5, vjust = 0.5, inherit.aes = FALSE, show.legend = FALSE) +

  annotate("text",
    x = 9.5, y = 9.65,
    label = "above null baseline\n→",
    hjust = 1, vjust = 1,
    size = 0.30 * FONT_SIZE, color = "#2166AC", fontface = "italic"
  ) +
  annotate("text",
    x = -4.8, y = 9.65,
    label = "below null baseline\n←",
    hjust = 0, vjust = 1,
    size = 0.30 * FONT_SIZE, color = "#B2182B", fontface = "italic"
  ) +

  annotate("segment",
    y = 5.5, yend = 5.5,
    x = -5, xend = 10,
    linewidth = 0.5, linetype = "dashed", color = "grey55"
  ) +
  annotate("text",
    x = 9.5, y = 5.7,
    label = "full-set null",
    hjust = 1, vjust = 0,
    size = 0.28 * FONT_SIZE, color = "black", fontface = "italic"
  ) +
  annotate("text",
    x = -4.8, y = 5.3,
    label = "5,000-array\nsubsample null",
    hjust = 0, vjust = 1,
    size = 0.28 * FONT_SIZE, color = "black", fontface = "italic"
  ) +

  scale_fill_manual(values = GROUP_COLORS, name = NULL) +
  coord_cartesian(xlim = c(-5, 10)) +
  scale_x_continuous(breaks = seq(-5, 10, by = 5), expand = expansion(mult = c(0.02, 0.02))) +

  labs(
    x = expression(Delta * " log"[2] * " (pseudocount-adjusted) CUT&Tag signal vs null baseline"),
    y = NULL,
    title = paste0("Relative repeat-level ", GROUP_LABEL, " signal")
  ) +

  theme_minimal(base_size = FONT_SIZE) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    axis.text.y = element_text(size = FONT_SIZE - 1, color = "black"),
    axis.text.x = element_text(size = FONT_SIZE - 1, color = "black"),
    axis.title.x = element_text(size = FONT_SIZE - 1, margin = margin(t = 8), color = "black"),
    plot.title = element_text(size = FONT_SIZE, face = "bold", margin = margin(b = 4)),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = FONT_SIZE - 1),
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(b = 4),
    plot.margin = margin(t = 8, r = 14, b = 6, l = 8)
  )

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir <- file.path(PERIOD_DIR, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".pdf")), p,
       width = 6.5, height = 9, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".png")), p,
       width = 6.5, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".svg")), p,
       width = 6.5, height = 9, dpi = 300, bg = "white", device = svglite::svglite)

cat("\nMedian Δ per bin (merged arrays):\n")
cat(sprintf("%-6s %-30s %8s %10s %10s\n", "Bin", "Label", "n_arrays", "rep1 Δ", "rep2 Δ"))
cat(strrep("-", 68), "\n")
for (b in ALL_BINS) {
  d1 <- bin_summary[sample == SAMPLES[1] & bin_id == b]
  d2 <- bin_summary[sample == SAMPLES[2] & bin_id == b]
  if (nrow(d1) > 0 && nrow(d2) > 0) {
    cat(sprintf("%-6s %-30s %8s %+10.3f %+10.3f\n",
        b, BIN_LABELS[as.character(b)], format(d1$n, big.mark = ","),
        d1$med_delta, d2$med_delta))
  }
}

message(sprintf("\nSaved: %s.pdf / .png / .svg", OUT_PREFIX))
