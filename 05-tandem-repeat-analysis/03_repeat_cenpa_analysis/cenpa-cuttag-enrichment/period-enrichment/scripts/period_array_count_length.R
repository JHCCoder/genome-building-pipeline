#!/usr/bin/env Rscript
# ============================================================================
# period_array_count_length.R — companion horizontal BARPLOTS for the merged-
#                    array period violin (period_violin_array.R): per TRF period
#                    bin, (a) the NUMBER of merged repeat arrays and (b) their
#                    MEAN LENGTH.
#
#   One horizontal bar per bin; the 9 bins are stacked in the SAME order as the
#   violin (bin 1, 1-10 bp microsatellites, at the BOTTOM; bin 9, 391+ bp, at
#   the TOP) so the three panels line up row-for-row when the strip is placed
#   to the LEFT of the violin.
#
#   Units = MERGED repeat arrays (bedtools merge -d 0 within each period bin) —
#   the same sampling unit as period_violin_array.R. Counts and lengths are
#   sample-independent (CUT&Tag signal is quantified at the same arrays).
#
#   Both axes are log10 (counts span 233-721,950; mean lengths span 60 bp to
#   142 kb, so a linear axis would squash the small bins into slivers). Value
#   labels are drawn at each bar end in human-readable form. Bar colors use the
#   same blues as the CENP-A enrichment violin (period_violin_array.R): dark
#   #1B4F72 = rep1 for the count panel, light #7FB3D8 = rep2 for the length
#   panel. X-axis titles read Log10(# of repeats) and Log10(mean repeat length).
#   Panel layout, font sizes and gridlines deliberately mirror
#   period_violin_array.R (FONT_SIZE knob, theme_minimal, x-gridlines only,
#   9 in tall) so the strip and violin assemble into one coherent figure.
#
#   The per-bin stats are computed from
#     ../data/merged/arrays/merged_arrays_all_bins.bed
#   and cached to ../results/period_array_count_length.tsv (reused when newer
#   than the BED), so re-rendering after a plot-format tweak is instant.
#
# Usage:
#   conda activate r-visualizations
#   Rscript period_array_count_length.R                    # default OUT_WIDTH x OUT_HEIGHT
#   Rscript period_array_count_length.R 8 9                # custom WIDTH HEIGHT (inches)
#   # omit y-axis bin labels + y-axis title (cleanest when assembled next to
#   # the violin, which already carries the bin names):
#   Rscript period_array_count_length.R --no-ylabels
#
# Output size: set OUT_WIDTH / OUT_HEIGHT (inches) below, or pass WIDTH HEIGHT
# as positional args (the notebook cell does this). Default 7.5 x 9 in (7.5
# gives the x-axis tick labels room not to overlap). OUT_HEIGHT 9 keeps the 9
# bin rows aligned with period_violin_array.R (both 9 in tall); changing
# OUT_HEIGHT means re-tuning TOP_MARGIN_PT / BOTTOM_MARGIN_PT below.
#
# Saved: ../plots/period_array_count_length.{pdf,png,svg}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

PERIOD_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
MERGED_BED <- file.path(PERIOD_DIR, "data", "merged", "arrays", "merged_arrays_all_bins.bed")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")
CACHE <- file.path(RESULTS_DIR, "period_array_count_length.tsv")

# ── Optional: --no-ylabels hides the bin names (for assembly next to violin) ──
args <- commandArgs(trailingOnly = TRUE)
SHOW_YLABELS <- !("--no-ylabels" %in% args)

# ── Font size (pt) — single knob, matches period_violin_array.R ───────────────
FONT_SIZE <- 12

# ── Output size (inches) — set to whatever you want to save out ───────────────
# Default 7.5 x 9 gives the x-axis tick labels (100, 1k, ...) enough room not
# to overlap at FONT_SIZE 12. OUT_HEIGHT 9 keeps the bin rows aligned with
# period_violin_array.R (changing it means re-tuning TOP_MARGIN_PT below).
# Two optional positional args override these, so the notebook cell can save at
# a custom size:  Rscript period_array_count_length.R WIDTH HEIGHT
OUT_WIDTH  <- 7.5
OUT_HEIGHT <- 9
pos_args <- args[!args %in% "--no-ylabels"]
if (length(pos_args) >= 1) OUT_WIDTH  <- as.numeric(pos_args[1])
if (length(pos_args) >= 2) OUT_HEIGHT <- as.numeric(pos_args[2])

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
ALL_BINS <- 1:9

# Bar colors match the CENP-A enrichment violin (period_violin_array.R):
# dark blue #1B4F72 = CENP-A rep1 -> count panel; light blue #7FB3D8 = rep2 -> length panel.
BAR_FILL_COUNT  <- "#1B4F72"
BAR_FILL_LENGTH <- "#7FB3D8"

# Top margin (pt) tuned so the 9 bin rows match the period_violin_array.R
# figure row-for-row when the strip is placed to its left (both figures are
# 9 in = 648 pt tall; same discrete y-axis expansion; strip bin pitch must
# equal the violin's ~59.07 pt -> panel height ~543 pt). The empty top band
# lines up with the violin's title/legend.
TOP_MARGIN_PT <- 51
BOTTOM_MARGIN_PT <- 4.1

# ── Data: per-bin count + mean length of merged repeat arrays ─────────────────
if (!file.exists(CACHE) || file.mtime(CACHE) < file.mtime(MERGED_BED)) {
  message("Computing per-bin array stats from ", MERGED_BED, " ...")
  bed <- fread(MERGED_BED, header = FALSE,
               col.names = c("chrom", "start", "end", "array_id", "bin_id", "n_intervals"))
  stats <- bed[, .(n_arrays = .N, mean_length = mean(end - start)), by = bin_id]
  stats[, bin_id := as.integer(bin_id)]
  setorder(stats, bin_id)
  fwrite(stats, CACHE, sep = "\t")
} else {
  stats <- fread(CACHE, sep = "\t")
}
stats <- stats[bin_id %in% ALL_BINS]
stats[, bin_label := factor(BIN_LABELS[as.character(bin_id)],
                            levels = BIN_LABELS[as.character(ALL_BINS)])]

# ── Value-label formatting ────────────────────────────────────────────────────
count_lab <- function(x) vapply(x, function(v) {
  if (v >= 1e6) sprintf("%.1fM", v / 1e6)          # 1.2M
  else if (v >= 1e5) sprintf("%.0fk", v / 1e3)      # 722k, 344k
  else if (v >= 1e3) sprintf("%.1fk", v / 1e3)      # 29.7k, 12.6k
  else sprintf("%.0f", v)                          # 799, 233
}, character(1))
bp_lab <- function(x) vapply(x, function(v) {
  if (v >= 100e3) sprintf("%.0f kb", v / 1e3)      # 142 kb
  else if (v >= 1e3) sprintf("%.1f kb", v / 1e3)    # 1.5 kb, 2.4 kb
  else sprintf("%.0f bp", v)                        # 60 bp, 838 bp
}, character(1))

# Bars drawn on a log10 axis floored at a round number below each bin's minimum
# (100 arrays; 10 bp), so bars start at the axis origin instead of log10(0).
FLOOR_C <- 2   # log10(100)  — count bars start at 100 arrays
FLOOR_L <- 1   # log10(10)   — length bars start at 10 bp
stats[, xc_count  := log10(n_arrays)   - FLOOR_C]
stats[, xc_length := log10(mean_length) - FLOOR_L]

# ── Panel 1: number of merged repeat arrays ───────────────────────────────────
p_count <- ggplot(stats, aes(y = bin_label, x = xc_count)) +
  geom_col(fill = BAR_FILL_COUNT, width = 0.68) +
  geom_text(aes(x = xc_count + 0.06, label = count_lab(n_arrays)),
            hjust = 0, size = 0.30 * FONT_SIZE, color = "black") +
  scale_x_continuous(
    breaks = 0:4,
    labels = c("100", "1k", "10k", "100k", "1M"),
    limits = c(0, 4.7),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  labs(x = expression("Log"[10] * "(# of repeats)"),
       y = if (SHOW_YLABELS) "Repeat array bin by period size" else NULL) +
  theme_minimal(base_size = FONT_SIZE) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    axis.text.y = if (SHOW_YLABELS) {
      element_text(size = FONT_SIZE, color = "black")
    } else {
      element_blank()
    },
    axis.ticks.y = if (SHOW_YLABELS) element_line(color = "black") else element_blank(),
    axis.text.x = element_text(size = FONT_SIZE, color = "black"),
    axis.title.x = element_text(size = FONT_SIZE, margin = margin(t = 8), color = "black"),
    axis.title.y = element_text(size = FONT_SIZE, color = "black"),
    plot.margin = margin(t = TOP_MARGIN_PT, r = 14, b = BOTTOM_MARGIN_PT, l = 8)
  )

# ── Panel 2: mean length of merged repeat arrays (no y-axis label) ────────────
p_length <- ggplot(stats, aes(y = bin_label, x = xc_length)) +
  geom_col(fill = BAR_FILL_LENGTH, width = 0.68) +
  geom_text(aes(x = xc_length + 0.06, label = bp_lab(mean_length)),
            hjust = 0, size = 0.30 * FONT_SIZE, color = "black") +
  scale_x_continuous(
    breaks = 0:5,
    labels = c("10", "100", "1k", "10k", "100k", "1M"),
    limits = c(0, 5.2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  labs(x = expression("Log"[10] * "(mean repeat length)")) +
  theme_minimal(base_size = FONT_SIZE) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = FONT_SIZE, color = "black"),
    axis.title.x = element_text(size = FONT_SIZE, margin = margin(t = 8), color = "black"),
    axis.title.y = element_blank(),
    plot.margin = margin(t = TOP_MARGIN_PT, r = 14, b = BOTTOM_MARGIN_PT, l = 8)
  )

# ── Assemble strip (count | mean length), same 9-in height as the violin ──────
p <- p_count + p_length + plot_layout(ncol = 2, widths = c(1, 1))

# ── Save ──────────────────────────────────────────────────────────────────────
out_dir <- file.path(PERIOD_DIR, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
OUT_PREFIX <- "period_array_count_length"

ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".pdf")), p,
       width = OUT_WIDTH, height = OUT_HEIGHT, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".png")), p,
       width = OUT_WIDTH, height = OUT_HEIGHT, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".svg")), p,
       width = OUT_WIDTH, height = OUT_HEIGHT, dpi = 300, bg = "white", device = svglite::svglite)

# ── Console summary ───────────────────────────────────────────────────────────
cat("\nMerged repeat arrays per period bin:\n")
cat(sprintf("%-6s %-30s %10s %12s\n", "Bin", "Label", "n_arrays", "mean_len_bp"))
cat(strrep("-", 60), "\n")
for (b in ALL_BINS) {
  r <- stats[bin_id == b]
  cat(sprintf("%-6s %-30s %10s %12.1f\n", b, BIN_LABELS[as.character(b)],
              format(r$n_arrays, big.mark = ","), r$mean_length))
}
message(sprintf("\nSaved: %s.pdf / .png / .svg", OUT_PREFIX))
