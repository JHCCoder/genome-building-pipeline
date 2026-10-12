#!/usr/bin/env Rscript
# ============================================================================
# period_violin_array.R — PLOTTING half of the horizontal per-merged-ARRAY
#                    Δ signal violin (Figure 4C for CENP-A): observed log2
#                    signal centered on the permutation null baseline, by TRF
#                    period bin.
#
#   DESCRIPTIVE layer. Each array's centered value is
#     Δ_i = log2(CPM_i + pseudo) − null_median(bin, sample)
#   where null_median = mean of the 1,000 per-iteration shuffled medians
#   (chr- and length-matched), computed on MERGED REPEAT ARRAYS (bedtools
#   merge -d 0 per bin) — the same unit as the observed values:
#     Bins 1–5: 1,000 iterations on a 5,000-array subsample
#     Bins 6–9: 1,000 iterations on all arrays (full set)
#   (dashed line between bins 5 and 6). White dot = median.
#   Significance stars are the FORMAL layer — BH-adjusted two-sided empirical
#   permutation P values (per sample) from the merged-array results CSV:
#   * < 0.05, ** < 0.01, *** < 0.001.
#
# Revised for Genome Research (October 2026): 8 pt Helvetica (Nimbus Sans),
# printed size 56 x 85 mm, no title, no y-axis bin labels (Figure 4B, placed to
# the left, carries them; both use the same 9-row discrete axis and the same
# plot-area top/bottom so rows line up), replicate key above the plot. The
# "above/below null baseline" arrows and the null-type notes moved to the
# figure legend.
#
# All expensive computation is done by scripts/period_violin_array_prepare.R,
# which caches:
#     $FIG4_PERIOD_DIR/results/period_violin_<group>_merged_delta.tsv
#     $FIG4_PERIOD_DIR/results/period_violin_<group>_merged_stars.tsv
# Star x/y positions come from the violin geometry (DODGE_W); if DODGE_W
# changes, re-run the _prepare script. FONT_SIZE / size changes never need it.
#
# Usage (r-visualizations env):
#   Rscript period_violin_array.R [CENP-A|H3K27ac] [WIDTH_MM HEIGHT_MM]
#   CENP-A  -> $FIG4_OUT_DIR/figure4C_cenpa_period_violin.{pdf,svg,png}
#   H3K27ac -> $FIG4_PERIOD_DIR/plots/period_violin_h3k27ac_merged.{pdf,svg,png}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

PERIOD_DIR  <- .cfg("FIG4_PERIOD_DIR")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

args <- commandArgs(trailingOnly = TRUE)
GROUP <- if (length(args) >= 1) args[1] else "CENP-A"
if (!GROUP %in% c("CENP-A", "H3K27ac")) stop("Unknown GROUP: ", GROUP)
FONT_SIZE <- 8                                       # pt at final print size; Genome Research: 8-10 pt
PANEL_W   <- (if (length(args) >= 2) as.numeric(args[2]) else 56) / 25.4   # in
PANEL_H   <- (if (length(args) >= 3) as.numeric(args[3]) else 85) / 25.4   # in
# Plot-area top/bottom matched to Figure 4A/4B (measured at 300 dpi)
TOP_MARGIN_PT    <- 7.1
BOTTOM_MARGIN_PT <- 7.0

GROUP_CONFIG <- list(
  "CENP-A"  = list(samples = c("XG_150", "XG_151"), label = "CENP-A",
                   colors = c("#1B4F72", "#7FB3D8"), prefix = "period_violin_cenpa_merged",
                   out = file.path(.cfg("FIG4_OUT_DIR"), "figure4C_cenpa_period_violin")),
  "H3K27ac" = list(samples = c("XG_152", "XG_153"), label = "H3K27ac",
                   colors = c("#B2182B", "#F4A582"), prefix = "period_violin_h3k27ac_merged",
                   out = file.path(PERIOD_DIR, "plots", "period_violin_h3k27ac_merged"))
)
cfg          <- GROUP_CONFIG[[GROUP]]
SAMPLES      <- cfg$samples
GROUP_LABEL  <- cfg$label
REP_LABELS   <- paste0(GROUP_LABEL, " rep ", 1:2)
GROUP_COLORS <- setNames(cfg$colors, REP_LABELS)

CONTROL_BINS <- 1:5
ALL_BINS     <- 1:9
DODGE_W      <- 0.8                                  # must match the _prepare script (star positions)

DELTA_FILE <- file.path(RESULTS_DIR, paste0(cfg$prefix, "_delta.tsv"))
STARS_FILE <- file.path(RESULTS_DIR, paste0(cfg$prefix, "_stars.tsv"))
if (!file.exists(DELTA_FILE)) {
  stop("Cached Δ-signal table not found: ", DELTA_FILE,
       "\n  Run scripts/period_violin_array_prepare.R ", GROUP, " first.")
}
dt <- fread(DELTA_FILE, sep = "\t")
sig_stars <- fread(STARS_FILE, sep = "\t")

dt[, bin_label := factor(as.character(bin_id), levels = as.character(ALL_BINS))]
dt[, sample_label := factor(REP_LABELS[match(sample, SAMPLES)], levels = REP_LABELS)]
bin_summary <- dt[, .(med_delta = median(delta, na.rm = TRUE), n = .N), by = .(sample, bin_id)]

p <- ggplot(dt, aes(x = delta, y = bin_label)) +
  geom_vline(xintercept = 0, linewidth = 0.5 / .pt * 1.5, color = "grey45") +
  geom_violin(aes(fill = sample_label), color = NA, alpha = 0.7, scale = "width", trim = TRUE,
              position = position_dodge(width = DODGE_W)) +
  stat_summary(aes(group = sample_label), fun = median, geom = "point", shape = 21, size = 1.1,
               fill = "white", color = "grey20", stroke = 0.3,
               position = position_dodge(width = DODGE_W)) +
  geom_text(data = sig_stars, aes(x = x, y = y, label = label), size = FONT_SIZE / .pt,
            family = "Nimbus Sans", color = "black", hjust = 0.5, vjust = 0.75,
            inherit.aes = FALSE, show.legend = FALSE) +
  # full-set null above, 5,000-array subsample null below (explained in the legend)
  annotate("segment", y = max(CONTROL_BINS) + 0.5, yend = max(CONTROL_BINS) + 0.5, x = -5, xend = 10,
           linewidth = 0.5 / .pt * 1.5, linetype = "dashed", color = "grey55") +
  scale_fill_manual(values = GROUP_COLORS, name = NULL) +
  coord_cartesian(xlim = c(-5, 10)) +
  scale_x_continuous(breaks = seq(-5, 10, by = 5), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = paste0("Δlog2 ", GROUP_LABEL, " vs null"), y = NULL) +
  theme_minimal(base_size = FONT_SIZE, base_family = "Nimbus Sans") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = FONT_SIZE, color = "black"),
    axis.title.x = element_text(size = FONT_SIZE, margin = margin(t = 2), color = "black"),
    axis.line.x = element_line(linewidth = 0.5 / .pt * 1.5, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.5 / .pt * 1.5, color = "black"),
    axis.ticks.length.x = unit(1.5, "pt"),
    legend.position = "top", legend.justification = "left", legend.direction = "horizontal",
    legend.text = element_text(size = FONT_SIZE),
    legend.key.height = unit(2.6, "mm"), legend.key.width = unit(3.5, "mm"),
    legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -2, -4),
    plot.margin = margin(t = TOP_MARGIN_PT, r = 3, b = BOTTOM_MARGIN_PT, l = 2))

nimbus <- "/usr/share/fonts/urw-base35/NimbusSans-"
nimbus_as_helvetica <- list(
  plain = list(alias = "Helvetica", file = paste0(nimbus, "Regular.otf")),
  bold = list(alias = "Helvetica", file = paste0(nimbus, "Bold.otf")),
  italic = list(alias = "Helvetica", file = paste0(nimbus, "Italic.otf")),
  bolditalic = list(alias = "Helvetica", file = paste0(nimbus, "BoldItalic.otf")))
svg_user_fonts <- list("sans" = nimbus_as_helvetica, "Nimbus Sans" = nimbus_as_helvetica,
                       "Helvetica" = nimbus_as_helvetica)

out_stem <- cfg$out
dir.create(dirname(out_stem), recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
for (ext in c("pdf", "svg", "png")) {
  f <- paste0(out_stem, ".", ext)
  if (file.exists(f)) file.rename(f, paste0(f, ".pre-", stamp))
}
ggsave(paste0(out_stem, ".pdf"), p, width = PANEL_W, height = PANEL_H, units = "in", device = cairo_pdf)
ggsave(paste0(out_stem, ".svg"), p, width = PANEL_W, height = PANEL_H, units = "in", user_fonts = svg_user_fonts)
ggsave(paste0(out_stem, ".png"), p, width = PANEL_W, height = PANEL_H, units = "in", dpi = 600, bg = "white")

cat("\nMedian Δ per bin (merged arrays):\n")
print(dcast(bin_summary, bin_id + n ~ sample, value.var = "med_delta")[order(bin_id)], digits = 3)
cat("Stars drawn:\n"); print(sig_stars[, .(bin_id, sample, label)])
cat(sprintf("Saved: %s.{pdf,svg,png} (%.0f x %.0f mm)\n", out_stem, PANEL_W * 25.4, PANEL_H * 25.4))
