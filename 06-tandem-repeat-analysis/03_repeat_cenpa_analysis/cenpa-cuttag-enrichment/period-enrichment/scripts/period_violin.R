#!/usr/bin/env Rscript
# ============================================================================
# period_violin.R — Horizontal violin: per-interval Δ signal (observed log2
#                    signal centered on the permutation null baseline) by TRF
#                    period bin, for a sample group.
#
#   DESCRIPTIVE layer. Each interval's centered value is
#     Δ_i = log2(CPM_i + pseudo) − null_median(bin, sample)
#   where null_median = mean of the 1,000 per-iteration shuffled medians
#   (chr- and length-matched):
#     Bins 1–5: 1,000 iterations on a 5,000-interval subsample
#     Bins 6–9: 1,000 iterations on all intervals (full set)
#   The violin uses ALL foreground intervals in each bin; the null already
#   marginalizes over fresh 5,000-interval subsample draws per iteration, so
#   centering on null_median is balanced for bins 1–5.
#   Δ is NOT an exact raw-CPM fold change: the pseudocount and the log-space
#   baseline make Δ a pseudocount-adjusted enrichment score relative to the
#   null baseline (2^B), not a literal 2×/4× fold change in CPM.
#   Significance stars are the FORMAL layer — empirical permutation p-values
#   (matched paired test for bins 1–5) read from the results CSV; they are not
#   derived from the per-interval values plotted here.
#
# Usage:
#   conda activate r-visualizations
#   Rscript period_violin.R CENP-A    # default; writes period_violin_cenpa.{pdf,png,svg}
#   Rscript period_violin.R H3K27ac   # writes period_violin_h3k27ac.{pdf,png,svg}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

PERIOD_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PARENT_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"

# ── Library sizes ─────────────────────────────────────────────────────────────
lib_file <- file.path(PARENT_DIR, "data", "qc", "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
} else {
  lib_sizes <- c("XG_150" = 27156788, "XG_151" = 25305681,
                 "XG_152" = 30805412, "XG_153" = 111250740)
}

# ── Sample group (CENP-A default; H3K27ac for the H3K27ac violin) ────────────
args <- commandArgs(trailingOnly = TRUE)
GROUP <- if (length(args) >= 1) args[1] else "CENP-A"
if (!GROUP %in% c("CENP-A", "H3K27ac")) stop("Unknown GROUP: ", GROUP)

GROUP_CONFIG <- list(
  "CENP-A"  = list(samples = c("XG_150", "XG_151"), label = "CENP-A",
                   colors = c("#1B4F72", "#7FB3D8"), prefix = "period_violin_cenpa"),
  "H3K27ac" = list(samples = c("XG_152", "XG_153"), label = "H3K27ac",
                   colors = c("#B2182B", "#F4A582"), prefix = "period_violin_h3k27ac")
)
cfg         <- GROUP_CONFIG[[GROUP]]
SAMPLES     <- cfg$samples
GROUP_LABEL <- cfg$label
GROUP_COLORS <- setNames(cfg$colors, paste0(GROUP_LABEL, "  rep", 1:2))
OUT_PREFIX  <- cfg$prefix

CENTROMERIC_BINS  <- 6:9
CONTROL_BINS      <- 1:5
ALL_BINS          <- c(CONTROL_BINS, CENTROMERIC_BINS)
DODGE_W           <- 0.8   # side-by-side dodge width for the two replicates

# ── Labels ───────────────────────────────────────────────────────────────────
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

# ── Load foreground ───────────────────────────────────────────────────────────
signal_list <- list()
for (s in SAMPLES) {
  f <- file.path(PERIOD_DIR, "data", "counts", paste0("trf_signal_", s, ".tsv"))
  x <- fread(f, header = FALSE,
    col.names = c("chrom", "start", "end", "interval_id",
                  "bin_id", "period_size", "copies_aligned",
                  "match_percent", "mean_coverage"))
  x[, sample := s]
  signal_list[[s]] <- x
}
signal <- rbindlist(signal_list)

# Normalize
signal[, lib_size := lib_sizes[sample]]
signal[, norm_signal := mean_coverage / (lib_size / 1e6)]
pseudo <- min(signal$norm_signal[signal$norm_signal > 0], na.rm = TRUE) / 2
signal[, log2_signal := log2(norm_signal + pseudo)]

# ── Compute null per bin from permutation background ─────────────────────────

# Control bins (1-5): fresh 5,000-subsample per iteration × 1,000 iterations
# (7 cols WITH header: bin_id, interval_id, iter, chrom, start, end, mean_signal).
# Null = mean of per-iteration medians (same definition as bins 6-9).
compute_null_ctrl <- function(bg_file, sample_name, bins) {
  if (!file.exists(bg_file)) return(NULL)
  x <- fread(bg_file, header = TRUE)
  setnames(x, "mean_signal", "mean_coverage")
  x[, sample := sample_name]
  x[, norm_signal := mean_coverage / (lib_sizes[sample_name] / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
  x <- x[bin_id %in% bins]
  bg_iter <- x[, .(bg_median = median(log2_signal, na.rm = TRUE)), by = .(bin_id, iter)]
  bg_iter[, .(null_median = mean(bg_median), n_iters = .N), by = bin_id]
}

# Full-set bins (6-9): all intervals × 1,000 iterations → null = mean of per-iteration medians
compute_null_cent <- function(bg_file, sample_name, bins) {
  if (!file.exists(bg_file)) return(NULL)
  x <- fread(bg_file, header = FALSE,
    col.names = c("chrom", "start", "end", "bin_id", "iter", "mean_coverage"))
  x[, sample := sample_name]
  x[, norm_signal := mean_coverage / (lib_sizes[sample_name] / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
  bg_iter <- x[bin_id %in% bins,
               .(bg_median = median(log2_signal, na.rm = TRUE)),
               by = .(bin_id, iter)]
  bg_iter[, .(null_median = mean(bg_median), n_iters = .N), by = bin_id]
}

null_all <- list()
for (s in SAMPLES) {
  ctrl_file <- file.path(PERIOD_DIR, "data", "counts",
                         paste0("trf_ctrl_bg_signal_", s, ".tsv"))
  null_ctrl <- compute_null_ctrl(ctrl_file, s, CONTROL_BINS)

  cent_file <- file.path(PERIOD_DIR, "data", "permutation", "counts",
                         paste0("trf_bg_signal_", s, ".tsv"))
  null_cent <- compute_null_cent(cent_file, s, CENTROMERIC_BINS)

  if (!is.null(null_ctrl)) {
    null_ctrl[, sample := s]
    null_all[[paste0(s, "_ctrl")]] <- null_ctrl
  }
  if (!is.null(null_cent)) {
    null_cent[, sample := s]
    null_all[[paste0(s, "_cent")]] <- null_cent
  }
}
null_dt <- rbindlist(null_all)

# ── Compute delta ─────────────────────────────────────────────────────────────
# Foreground uses ALL intervals in each bin (null for bins 1-5 already
# marginalizes over fresh 5,000-subsample draws per iteration).
dt <- signal[bin_id %in% ALL_BINS]
dt <- merge(dt, null_dt[, .(sample, bin_id, null_median)], by = c("sample", "bin_id"))
dt[, delta := log2_signal - null_median]

# Bin sizes: per-bin interval counts (identical across samples — both measure
# signal at the same TRF intervals; interval_id is globally unique).
bin_n <- signal[bin_id %in% ALL_BINS, .(n = uniqueN(interval_id)), by = bin_id]
BIN_LABELS_N <- setNames(
  vapply(ALL_BINS, function(b) {
    paste0(BIN_LABELS[as.character(b)],
           "\n(n=", format(bin_n[bin_id == b, n], big.mark = ","), ")")
  }, character(1)),
  as.character(ALL_BINS)
)

# Bin labels (reversed — bin 1 at top)
dt[, bin_label := factor(BIN_LABELS_N[as.character(bin_id)],
                          levels = rev(BIN_LABELS_N[as.character(ALL_BINS)]))]

# Sample labels
dt[, sample_label := paste0(GROUP_LABEL, "  rep", match(sample, SAMPLES))]

# Per-bin summary
bin_summary <- dt[, .(med_delta = median(delta, na.rm = TRUE), n = .N),
                  by = .(sample, bin_id, bin_label)]

# ── Significance stars ────────────────────────────────────────────────────────
# Per-replicate significance read from the permutation-test results CSV written
# by 3_analyze_period_enrichment.R. The `sig` column already encodes the star
# count (*, **, ***) so the figure shows the real number of asterisks.
# Placement: each star sits just right of ITS OWN violin's visible right tail
# (the drawn polygon's last point where width > 2% of that violin's max width,
# + a small offset), at the vertical centre of that replicate's dodged violin
# band (rep1 lower, rep2 upper) — so it clearly hugs its own violin and is
# unambiguous per replicate. Anchoring to the drawn tail (not a high percentile)
# keeps stars tight against the violin even when extreme outliers would inflate
# a percentile.
res_file <- file.path(PERIOD_DIR, "results", "period_enrichment_permutation_results.csv")
sig_stars <- data.table(bin_id = integer(0), sample = character(0))
if (file.exists(res_file)) {
  res <- fread(res_file)
  sig_stars <- res[sample %in% SAMPLES & p_adj < 0.05,
                   .(bin_id, sample, label = sig)]
}
if (nrow(sig_stars) > 0) {
  # Build the violin geometry once to read each violin's drawn right tail.
  probe <- ggplot(dt, aes(x = delta, y = bin_label, fill = sample_label)) +
    geom_violin(color = NA, alpha = 0.55, scale = "width", trim = TRUE,
                position = position_dodge(width = DODGE_W)) +
    scale_fill_manual(values = GROUP_COLORS, name = NULL)
  vg <- as.data.table(ggplot_build(probe)$data[[1]])
  # Map each violin group to (sample, bin_id) via its y-centre.
  # rep1 (lower) sits at bin_pos - DODGE_W/4, rep2 (upper) at bin_pos + DODGE_W/4;
  # bin_pos = 10 - bin_id (bin 1 top, bin 9 bottom on the reversed factor axis).
  gmap <- vg[, .(yc = mean(y)), by = group]
  gmap[, bin_pos := round(yc)]
  gmap[, sample := ifelse(yc < bin_pos, SAMPLES[1], SAMPLES[2])]
  gmap[, bin_id := 10 - bin_pos]
  gmap[, tail_edge := vapply(group, function(g) {
    sub <- vg[group == g]
    mw <- max(sub$violinwidth)
    max(sub[violinwidth > 0.02 * mw, x])
  }, numeric(1))]
  sig_stars <- merge(sig_stars, gmap[, .(sample, bin_id, tail_edge)],
                     by = c("sample", "bin_id"))
  # y: centre of that replicate's dodged violin band (DODGE_W/4, verified with
  # ggplot_build); x: just right of the visible right tail.
  sig_stars[, y := (10 - bin_id) + ifelse(sample == SAMPLES[1], -DODGE_W/4, DODGE_W/4)]
  sig_stars[, x := tail_edge + 0.3]
} else {
  sig_stars[, `:=`(x = numeric(), y = numeric(), label = character())]
}

# ── Plot ─────────────────────────────────────────────────────────────────────
p <- ggplot(dt, aes(x = delta, y = bin_label)) +

  # Zero reference
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.45, color = "grey45") +

  # Violin (dodged so the two replicates sit side by side)
  geom_violin(
    aes(fill = sample_label),
    color = NA,
    alpha = 0.55,
    scale = "width",
    trim = TRUE,
    position = position_dodge(width = DODGE_W)
  ) +

  # Per-replicate median marker
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

  # Significance stars (just right of each violin's peak)
  geom_text(data = sig_stars, aes(x = x, y = y, label = label),
            size = 4.5, fontface = "bold", color = "black",
            hjust = 0.5, vjust = 0.5, inherit.aes = FALSE, show.legend = FALSE) +

  # Enrichment / depletion labels — tucked just above the top violin, inside the
  # discrete axis's natural top margin. (Previously placed at y = 10.7, which
  # pushed the y-axis top ~1.3 units away from the data and left a dead gap
  # between the legend and the "1–10 bp" row.) "above null" sits at the far
  # right where the top violin's tail is negligible; "below null" sits at the
  # far left where the top violin (min Δ ≈ +2) never reaches.
  annotate("text",
    x = 9.5, y = 9.45,
    label = "→  above null baseline",
    hjust = 1, vjust = 1,
    size = 3.2, color = "#2166AC", fontface = "italic"
  ) +
  annotate("text",
    x = -3.5, y = 9.45,
    label = "below null baseline  ←",
    hjust = 0, vjust = 1,
    size = 3.2, color = "#B2182B", fontface = "italic"
  ) +

  # Separator: subsample-null bins (1-5) above, full-set-null bins (6-9) below.
  # Bin 5 (196-347 bp) sits at y=5, bin 6 (348-349 bp) at y=4, so the divider
  # goes between them at y=4.5 (was wrongly at y=5.5 = between bins 4 and 5).
  annotate("segment",
    y = 4.5, yend = 4.5,
    x = -5, xend = 10,
    linewidth = 0.5, linetype = "dashed", color = "grey55"
  ) +
  annotate("text",
    x = 9.5, y = 4.35,
    label = "full-set null",
    hjust = 1, vjust = 1,
    size = 3, color = "black", fontface = "italic"
  ) +
  annotate("text",
    x = -4.8, y = 4.65,
    label = "subsample null",
    hjust = 0, vjust = 0,
    size = 3, color = "black", fontface = "italic"
  ) +

  # Scales
  scale_fill_manual(
    values = GROUP_COLORS,
    name = NULL
  ) +
  coord_cartesian(xlim = c(-5, 10)) +
  scale_x_continuous(
    breaks = seq(-5, 10, by = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  # Labels
  labs(
    x = expression(Delta * " log"[2] * " (pseudocount-adjusted) CUT&Tag signal vs null baseline"),
    y = NULL,
    title = paste0(GROUP_LABEL, " signal enrichment")
  ) +

  # Theme
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 9.5, margin = margin(t = 8), color = "black"),
    plot.title = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 9.5),
    legend.key.size = unit(0.6, "cm"),
    legend.margin = margin(b = 4),
    plot.margin = margin(t = 8, r = 14, b = 6, l = 8)
  )

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir <- file.path(PERIOD_DIR, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".pdf")), p,
       width = 6, height = 8, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".png")), p,
       width = 6, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(OUT_PREFIX, ".svg")), p,
       width = 6, height = 8, dpi = 300, bg = "white", device = svg)

# ── Print summary ────────────────────────────────────────────────────────────
cat("\nMedian Δ per bin:\n")
cat(sprintf("%-6s %-30s %8s %10s %10s\n", "Bin", "Label", "n", "rep1 Δ", "rep2 Δ"))
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
