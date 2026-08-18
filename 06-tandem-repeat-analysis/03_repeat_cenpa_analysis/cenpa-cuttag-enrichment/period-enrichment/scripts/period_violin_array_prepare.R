#!/usr/bin/env Rscript
# ============================================================================
# period_violin_array_prepare.R — DATA-PROCESSING half of the horizontal
#                   per-merged-ARRAY Δ signal violin (see period_violin_array.R
#                   for the figure description).
#
# This script contains ALL of the expensive computation that used to live
# inside period_violin_array.R:
#   * reads the merged-array signal files for the requested sample GROUP and
#     computes the per-array log2 signal + global pseudocount (group-specific)
#   * reads the ~1 GB ctrl (bins 1-5) and full (bins 6-9) permutation
#     background files per sample and computes the per-bin null_median
#     (mean of the 1,000 per-iteration shuffled medians)
#   * computes Δ_i = log2_signal − null_median(bin, sample) per merged array
#   * computes the significance-star positions (x = violin tail edge, y =
#     dodge-centred bin) via ggplot_build on a probe violin
#
# It writes two caches under ../results/, suffixed by the GROUP:
#     period_violin_<group>_merged_delta.tsv   (array_id, bin_id, sample, delta)
#     period_violin_<group>_merged_stars.tsv   (bin_id, sample, label, x, y)
#
# The PLOTTING half is scripts/period_violin_array.R, which reads these two
# tables and renders the figure in seconds. Run THIS script only when the raw
# data changes; for any pure plotting tweak just re-run the plot script.
#
# NOTE: star x/y positions are part of the cache because they are computed from
# the violin geometry. If you change DODGE_W or anything that shifts the
# violin/dodge layout in the plot script, re-run this _prepare script so the
# star positions stay aligned. FONT_SIZE changes never need a re-run.
#
# Usage (run rarely — only when raw data changes):
#   conda activate r-visualizations
#   Rscript period_violin_array_prepare.R CENP-A
#   Rscript period_violin_array_prepare.R H3K27ac
#
# Saved: ../results/period_violin_<group>_merged_{delta,stars}.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

PERIOD_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PARENT_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
MERGED_DIR <- file.path(PERIOD_DIR, "data", "merged")
COUNTS_DIR <- file.path(MERGED_DIR, "counts")
BG_COUNTS_DIR <- file.path(MERGED_DIR, "permutation", "counts")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

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
                   colors = c("#1B4F72", "#7FB3D8"), prefix = "period_violin_cenpa_merged"),
  "H3K27ac" = list(samples = c("XG_152", "XG_153"), label = "H3K27ac",
                   colors = c("#B2182B", "#F4A582"), prefix = "period_violin_h3k27ac_merged")
)
cfg         <- GROUP_CONFIG[[GROUP]]
SAMPLES     <- cfg$samples
GROUP_LABEL <- cfg$label
GROUP_COLORS <- setNames(cfg$colors, paste0(GROUP_LABEL, "  rep", 1:2))

CENTROMERIC_BINS  <- 6:9
CONTROL_BINS      <- 1:5
ALL_BINS          <- c(CONTROL_BINS, CENTROMERIC_BINS)
DODGE_W           <- 0.8

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

# ── Load foreground (merged arrays) ───────────────────────────────────────────
signal_list <- list()
for (s in SAMPLES) {
  f <- file.path(COUNTS_DIR, paste0("trf_signal_array_", s, ".tsv"))
  if (!file.exists(f)) stop("Missing merged-array signal file: ", f)
  x <- fread(f, header = FALSE,
    col.names = c("chrom", "start", "end", "array_id",
                  "bin_id", "n_intervals", "mean_coverage"))
  x[, sample := s]
  signal_list[[s]] <- x
}
signal <- rbindlist(signal_list)

# Normalize
signal[, lib_size := lib_sizes[sample]]
signal[, norm_signal := mean_coverage / (lib_size / 1e6)]
pseudo <- min(signal$norm_signal[signal$norm_signal > 0], na.rm = TRUE) / 2
signal[, log2_signal := log2(norm_signal + pseudo)]

# ── Compute null per bin from permutation background (merged arrays) ──────────

# Control bins (1-5): fresh 5,000-array subsample per iteration × 1,000 iters
# (7 cols WITH header: bin_id, array_id, iter, chrom, start, end, mean_signal).
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

# Full-set bins (6-9): all arrays × 1,000 iterations.
# (7 cols no header: chrom, start, end, array_id, bin_id, iter, mean_coverage)
compute_null_cent <- function(bg_file, sample_name, bins) {
  if (!file.exists(bg_file)) return(NULL)
  x <- fread(bg_file, header = FALSE,
    col.names = c("chrom", "start", "end", "array_id", "bin_id", "iter", "mean_coverage"))
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
  ctrl_file <- file.path(COUNTS_DIR, paste0("trf_ctrl_bg_signal_array_", s, ".tsv"))
  null_ctrl <- compute_null_ctrl(ctrl_file, s, CONTROL_BINS)

  cent_file <- file.path(BG_COUNTS_DIR, paste0("trf_bg_signal_array_", s, ".tsv"))
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

# ── Compute delta (per merged array) ──────────────────────────────────────────
dt <- signal[bin_id %in% ALL_BINS]
dt <- merge(dt, null_dt[, .(sample, bin_id, null_median)], by = c("sample", "bin_id"))
dt[, delta := log2_signal - null_median]

# Bin sizes: per-bin MERGED-ARRAY counts (identical across samples — signal is
# measured at the same arrays).
bin_n <- signal[bin_id %in% ALL_BINS, .(n = uniqueN(array_id)), by = bin_id]
BIN_LABELS_N <- setNames(
  vapply(ALL_BINS, function(b) {
    paste0(BIN_LABELS[as.character(b)],
           "\n(n=", format(bin_n[bin_id == b, n], big.mark = ","), " arrays)")
  }, character(1)),
  as.character(ALL_BINS)
)

dt[, bin_label := factor(BIN_LABELS_N[as.character(bin_id)],
                          levels = BIN_LABELS_N[as.character(ALL_BINS)])]
dt[, sample_label := paste0(GROUP_LABEL, "  rep", match(sample, SAMPLES))]

bin_summary <- dt[, .(med_delta = median(delta, na.rm = TRUE), n = .N),
                  by = .(sample, bin_id, bin_label)]

# ── Significance stars (positions from violin geometry; label from results) ───
res_file <- file.path(PERIOD_DIR, "results", "period_enrichment_merged_permutation_results.csv")
sig_stars <- data.table(bin_id = integer(0), sample = character(0))
if (file.exists(res_file)) {
  res <- fread(res_file)
  sig_stars <- res[sample %in% SAMPLES & p_adj < 0.05,
                   .(bin_id, sample, label = sig)]
}
if (nrow(sig_stars) > 0) {
  probe <- ggplot(dt, aes(x = delta, y = bin_label, fill = sample_label)) +
    geom_violin(color = NA, alpha = 0.55, scale = "width", trim = TRUE,
                position = position_dodge(width = DODGE_W)) +
    scale_fill_manual(values = GROUP_COLORS, name = NULL)
  vg <- as.data.table(ggplot_build(probe)$data[[1]])
  gmap <- vg[, .(yc = mean(y)), by = group]
  gmap[, bin_pos := round(yc)]
  gmap[, sample := ifelse(yc < bin_pos, SAMPLES[1], SAMPLES[2])]
  gmap[, bin_id := bin_pos]
  gmap[, tail_edge := vapply(group, function(g) {
    sub <- vg[group == g]
    mw <- max(sub$violinwidth)
    max(sub[violinwidth > 0.02 * mw, x])
  }, numeric(1))]
  sig_stars <- merge(sig_stars, gmap[, .(sample, bin_id, tail_edge)],
                     by = c("sample", "bin_id"))
  sig_stars[, y := bin_id + ifelse(sample == SAMPLES[1], -DODGE_W/4, DODGE_W/4)]
  sig_stars[, x := tail_edge + 0.3]
} else {
  sig_stars[, `:=`(x = numeric(), y = numeric(), label = character())]
}

# ── Print the same per-bin summary table as the combined script ───────────────
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

# ── Write the caches the plot script reads ────────────────────────────────────
# delta cache: per-array Δ signal (array_id kept so the plot script can
# recompute the per-bin array counts exactly as the analysis does).
fwrite(dt[, .(array_id, bin_id, sample, delta)], DELTA_FILE, sep = "\t")
# stars cache: pre-computed star positions + labels.
fwrite(sig_stars[, .(bin_id, sample, label, x, y)], STARS_FILE, sep = "\t")
message(sprintf("\nSaved caches: %s (%d rows), %s (%d rows)",
                DELTA_FILE, nrow(dt), STARS_FILE, nrow(sig_stars)))
