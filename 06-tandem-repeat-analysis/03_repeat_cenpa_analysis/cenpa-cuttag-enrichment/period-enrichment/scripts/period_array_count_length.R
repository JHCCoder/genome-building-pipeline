#!/usr/bin/env Rscript
# ============================================================================
# period_array_count_length.R — Figure 4B: per TRF period bin, (left) the NUMBER
#                    of merged repeat arrays and (right) their MEAN LENGTH.
#
#   One horizontal bar per bin; the 9 bins are stacked in the SAME order as the
#   CENP-A violin (period_violin_array.R, Figure 4C): bin 1 (1-10 bp) at the
#   BOTTOM, bin 9 (391+ bp) at the TOP, so B and C line up row-for-row when C is
#   placed to the right of B (same height and same top/bottom plot margins).
#
#   Units = MERGED repeat arrays (bedtools merge -d 0 within each period bin) —
#   the same sampling unit as period_violin_array.R.
#
#   Both axes are log10 (counts span 233-721,950; mean lengths 60 bp-142 kb);
#   value labels are drawn at each bar end. Revised for Genome Research
#   (October 2026): 8 pt Helvetica (Nimbus Sans), printed size 58 x 85 mm, the
#   four bins that recur in Figure 4 (1-10, 193-195, 348-349, 386-390 bp) in
#   their Figure 4 colours (config.sh FIG4_COLOR_*) and the rest grey.
#
#   The per-bin stats are computed from
#     $FIG4_PERIOD_DIR/data/merged/arrays/merged_arrays_all_bins.bed
#   and cached to $FIG4_PERIOD_DIR/results/period_array_count_length.tsv
#   (reused when newer than the BED).
#
# Usage (r-visualizations env):
#   Rscript period_array_count_length.R [WIDTH_MM HEIGHT_MM]
#
# Saved: $FIG4_OUT_DIR/figure4B_period_array_count_length.{pdf,svg,png}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

.cfg <- function(name) {
  v <- Sys.getenv(name)
  if (nzchar(v)) return(v)
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  root <- normalizePath(if (length(arg)) dirname(sub("^--file=", "", arg[1])) else getwd())
  # repo-root config.sh (next to .git): this tree also has a local cenpa-cuttag-enrichment/config.sh
  while (!(file.exists(file.path(root, "config.sh")) && file.exists(file.path(root, ".git"))) &&
         root != dirname(root)) root <- dirname(root)
  cmd <- sprintf("source \"%s\" 2>/dev/null; printf \"%%s\" \"${%s}\"",
                 file.path(root, "config.sh"), name)
  v <- system2("bash", c("-c", shQuote(cmd)), stdout = TRUE)
  if (length(v) == 0 || !nzchar(v)) stop(sprintf("ERROR: %s is not set in the environment or in config.sh", name))
  v
}

PERIOD_DIR <- .cfg("FIG4_PERIOD_DIR")
MERGED_BED <- file.path(PERIOD_DIR, "data", "merged", "arrays", "merged_arrays_all_bins.bed")
CACHE      <- file.path(PERIOD_DIR, "results", "period_array_count_length.tsv")

args <- commandArgs(trailingOnly = TRUE)
FONT_SIZE <- 8                                       # pt at final print size; Genome Research: 8-10 pt
PANEL_W   <- (if (length(args) >= 1) as.numeric(args[1]) else 58) / 25.4   # in
PANEL_H   <- (if (length(args) >= 2) as.numeric(args[2]) else 85) / 25.4   # in
# Plot-area margins shared with period_violin_array.R so the 9 rows align (pt)
TOP_MARGIN_PT    <- 15.8   # pt; puts the plot-area top level with Figure 4A (measured), so A, B, C share one row
BOTTOM_MARGIN_PT <- 2

BIN_LABELS <- c(
  "1" = "1–10 bp",
  "2" = "11–50 bp",
  "3" = "51–192 bp",
  "4" = "193–195 bp",
  "5" = "196–347 bp",
  "6" = "348–349 bp",
  "7" = "350–385 bp",
  "8" = "386–390 bp",
  "9" = "391+ bp"
)
ALL_BINS <- 1:9
GREY <- "#9A9A9A"
BIN_COLORS <- c("1" = .cfg("FIG4_COLOR_1_10"), "2" = GREY, "3" = GREY,
                "4" = .cfg("FIG4_COLOR_193_195"), "5" = GREY,
                "6" = .cfg("FIG4_COLOR_348_349"), "7" = GREY,
                "8" = .cfg("FIG4_COLOR_386_390"), "9" = GREY)

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
stats[, fill := BIN_COLORS[as.character(bin_id)]]

count_lab <- function(x) vapply(x, function(v) {
  if (v >= 1e5) sprintf("%.0fk", v / 1e3)
  else if (v >= 1e3) sprintf("%.1fk", v / 1e3)
  else sprintf("%.0f", v)
}, character(1))
bp_lab <- function(x) vapply(x, function(v) {
  if (v >= 100e3) sprintf("%.0f kb", v / 1e3)
  else if (v >= 1e3) sprintf("%.1f kb", v / 1e3)
  else sprintf("%.0f bp", v)
}, character(1))

# Bars start at a round floor (100 arrays; 10 bp) on the log10 axis
FLOOR_C <- 2
FLOOR_L <- 1
stats[, xc_count  := log10(n_arrays)    - FLOOR_C]
stats[, xc_length := log10(mean_length) - FLOOR_L]

base_theme <- theme_minimal(base_size = FONT_SIZE, base_family = "Nimbus Sans") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    axis.text.x = element_text(size = FONT_SIZE, color = "black"),
    axis.title.x = element_text(size = FONT_SIZE, margin = margin(t = 2), color = "black"),
    axis.line.x = element_line(linewidth = 0.5 / .pt * 1.5, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.5 / .pt * 1.5, color = "black"),
    axis.ticks.length.x = unit(1.5, "pt"))

bar_panel <- function(xcol, labfun, valcol, breaks, labels, lim, xlab, show_y) {
  ggplot(stats, aes(y = bin_label, x = .data[[xcol]])) +
    geom_col(aes(fill = fill), width = 0.68) +
    scale_fill_identity() +
    geom_text(aes(x = .data[[xcol]] + 0.08, label = labfun(.data[[valcol]])),
              hjust = 0, size = FONT_SIZE / .pt, family = "Nimbus Sans", color = "black") +
    scale_x_continuous(breaks = breaks, labels = labels, limits = c(0, lim),
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    labs(x = xlab, y = NULL) +
    base_theme +
    theme(axis.text.y = if (show_y) element_text(size = FONT_SIZE, color = "black") else element_blank(),
          plot.margin = margin(t = TOP_MARGIN_PT, r = 2, b = BOTTOM_MARGIN_PT, l = 1))
}

# Counts: ticks at 100, 10k, 1M; lengths: 10 bp, 1 kb, 100 kb (every other decade, so 8 pt labels fit)
p_count  <- bar_panel("xc_count", count_lab, "n_arrays", c(0, 2, 4), c("100", "10k", "1M"), 5.4,
                      "Merged arrays", TRUE)
p_length <- bar_panel("xc_length", bp_lab, "mean_length", c(0, 2, 4), c("10", "1k", "100k"), 6.4,
                      "Mean length (bp)", FALSE)
p <- p_count + p_length + plot_layout(ncol = 2, widths = c(1, 1))

# ── Save (back up existing outputs first) ─────────────────────────────────────
nimbus <- "/usr/share/fonts/urw-base35/NimbusSans-"
nimbus_as_helvetica <- list(
  plain = list(alias = "Helvetica", file = paste0(nimbus, "Regular.otf")),
  bold = list(alias = "Helvetica", file = paste0(nimbus, "Bold.otf")),
  italic = list(alias = "Helvetica", file = paste0(nimbus, "Italic.otf")),
  bolditalic = list(alias = "Helvetica", file = paste0(nimbus, "BoldItalic.otf")))
svg_user_fonts <- list("sans" = nimbus_as_helvetica, "Nimbus Sans" = nimbus_as_helvetica,
                       "Helvetica" = nimbus_as_helvetica)

out_dir <- .cfg("FIG4_OUT_DIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_stem <- file.path(out_dir, "figure4B_period_array_count_length")
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
for (ext in c("pdf", "svg", "png")) {
  f <- paste0(out_stem, ".", ext)
  if (file.exists(f)) file.rename(f, paste0(f, ".pre-", stamp))
}
ggsave(paste0(out_stem, ".pdf"), p, width = PANEL_W, height = PANEL_H, units = "in", device = cairo_pdf)
ggsave(paste0(out_stem, ".svg"), p, width = PANEL_W, height = PANEL_H, units = "in", user_fonts = svg_user_fonts)
ggsave(paste0(out_stem, ".png"), p, width = PANEL_W, height = PANEL_H, units = "in", dpi = 600, bg = "white")

cat("\nMerged repeat arrays per period bin:\n")
print(stats[, .(bin_id, bin = as.character(bin_label), n_arrays, mean_length = round(mean_length, 1))])
cat(sprintf("Saved: %s.{pdf,svg,png} (%.0f x %.0f mm)\n", out_stem, PANEL_W * 25.4, PANEL_H * 25.4))
