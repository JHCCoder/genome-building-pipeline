#!/usr/bin/env Rscript
# ============================================================================
# supp_monomer_bins_foreground_vs_null_overlay_array_prepare.R
# DATA-PROCESSING half of the merged-array per-unit overlay for the three CENP-A
# monomer families: bins 4 (193-195 bp), 6 (348-349 bp) and 8 (386-390 bp).
#
# This script contains ALL of the expensive computation that used to live
# inside supp_monomer_bins_foreground_vs_null_overlay_array.R:
#   * global pseudocount from all 4 merged samples (foreground signal files)
#   * foreground per-array log2 signal for bins 4/6/8, CENP-A samples
#   * null per-array log2 signal (pooled across 1,000 shuffled iterations),
#     stream-filtered from the ~1 GB background files; nulls >1M rows thinned
#     to 1M (set.seed 42) — the thinning is part of the cached data.
#   * descriptive observed vs null medians and panel direction/n_sig printouts
#     (direction read from the MERGED permutation results CSV).
#
# It writes the plotting data to
#     ../results/supp_monomer_bins_foreground_vs_null_distribution_merged_data.tsv
#     columns: bin_id, sample, grp ("foreground"/"null"), value (log2 normalized)
#
# The PLOTTING half is scripts/supp_monomer_bins_foreground_vs_null_overlay_array.R,
# which reads this table (plus the small merged permutation results CSV for the
# panel labels) and renders the figure in seconds. Run THIS script only when the
# raw data (signal/permutation counts, library sizes, binning) changes; for any
# pure plotting tweak just re-run the plot script.
#
# Normalization matches 3_analyze_period_enrichment_array.R EXACTLY (verified
# against results/period_enrichment_merged_permutation_results.csv):
#   norm_signal = mean_coverage / (lib_size / 1e6)
#   pseudo      = min(norm_signal > 0) / 2  (global, all 4 merged samples)
#   log2_signal = log2(norm_signal + pseudo)
#
# NOTE on the zero-coverage spike at log2 = -9.38: for a short bin (e.g. bin 4,
# XG_151) most shuffled arrays do not overlap CENP-A peaks, so >50% have mean
# coverage = 0, and log2(0 + pseudo) puts them on the spike at log2 = -9.38.
# This is a real property of the shuffled background (matches bg_mean_of_medians
# in the merged permutation CSV), not an artifact of the figure.
#
# Usage (run rarely — only when raw data changes):
#   conda activate r-visualizations
#   Rscript supp_monomer_bins_foreground_vs_null_overlay_array_prepare.R
#
# Saved: ../results/supp_monomer_bins_foreground_vs_null_distribution_merged_data.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

PERIOD_DIR  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PARENT_DIR  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
MERGED_DIR  <- file.path(PERIOD_DIR, "data", "merged")
COUNTS_DIR  <- file.path(MERGED_DIR, "counts")
BG_COUNTS_DIR <- file.path(MERGED_DIR, "permutation", "counts")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

SAMPLES <- c("XG_150", "XG_151")                       # CENP-A replicates
ALL_SAMPLES   <- c("XG_150", "XG_151", "XG_152", "XG_153")
BINS <- c(4, 6, 8)   # 193-195 bp, 348-349 bp, 386-390 bp

DATA_FILE <- file.path(RESULTS_DIR,
                       "supp_monomer_bins_foreground_vs_null_distribution_merged_data.tsv")

# ── Library sizes (must match the analysis) ────────────────────────────────
lib_file <- file.path(PARENT_DIR, "data", "qc", "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
} else {
  stop("Library size file not found: ", lib_file)
}

# ── Global pseudocount (all 4 merged samples — same as the analysis) ────────
sig_all <- rbindlist(lapply(ALL_SAMPLES, function(s) {
  f <- file.path(COUNTS_DIR, paste0("trf_signal_array_", s, ".tsv"))
  if (!file.exists(f)) stop("Missing merged-array signal file: ", f)
  x <- fread(f, header = FALSE,
             col.names = c("chrom", "start", "end", "array_id", "bin_id",
                           "n_intervals", "mean_coverage"))
  x[, sample := s]
}))
sig_all[, norm_signal := mean_coverage / (lib_sizes[sample] / 1e6)]
pseudo <- min(sig_all$norm_signal[sig_all$norm_signal > 0], na.rm = TRUE) / 2
message(sprintf("Global pseudocount (merged arrays): %.6e (log2 = %.3f)", pseudo, log2(pseudo)))

lgsig <- function(cov, s) log2(cov / (lib_sizes[s] / 1e6) + pseudo)

# ── Foreground: per-array signal for bins 4/6/8, CENP-A samples ─────────────
fg_all <- sig_all[bin_id %in% BINS & sample %in% SAMPLES]
fg_all[, value := lgsig(mean_coverage, sample)]

# ── Null per-array (pooled across 1,000 iterations) ─────────────────────────
# bin 4 -> ctrl null (bins 1-5: 5,000-array subsample/iter), 7 cols WITH header.
#   awk keeps the header line + bin_id == 4 rows (streamed, never fully loaded).
# bins 6 & 8 -> full-set null (bins 6-9), 7 cols NO header; awk keeps 6/8.
load_null <- function(s, b) {
  if (b %in% 1:5) {
    f <- file.path(COUNTS_DIR, paste0("trf_ctrl_bg_signal_array_", s, ".tsv"))
    cmd <- sprintf("awk -F'\\t' 'NR==1 || $1==%d' '%s'", b, f)
    x <- fread(cmd = cmd, header = TRUE)
    setnames(x, "mean_signal", "mean_coverage")
  } else {
    f <- file.path(BG_COUNTS_DIR, paste0("trf_bg_signal_array_", s, ".tsv"))
    cmd <- sprintf("awk -F'\\t' '$5==%d' '%s'", b, f)
    x <- fread(cmd = cmd, header = FALSE,
               col.names = c("chrom", "start", "end", "array_id", "bin_id", "iter", "mean_coverage"))
  }
  x
}

d_list <- list()
for (b in BINS) {
  for (s in SAMPLES) {
    fg <- fg_all[bin_id == b & sample == s]
    nu <- load_null(s, b)
    if (nrow(nu) > 1e6) { set.seed(42); nu <- nu[sample(.N, 1e6)] }  # thin ONLY for plotting
    d_list[[paste(b, s, "fg")]] <- data.table(bin_id = b, sample = s, grp = "foreground",
                                              value = lgsig(fg$mean_coverage, s))
    d_list[[paste(b, s, "nu")]] <- data.table(bin_id = b, sample = s, grp = "null",
                                              value = lgsig(nu$mean_coverage, s))
  }
}
d <- rbindlist(d_list)

# ── Descriptive medians + direction (same printout as the combined script) ──
SAMPLE_LABELS <- c("XG_150" = "CENP-A rep1", "XG_151" = "CENP-A rep2")
BIN_LABELS <- c("4" = "193-195 bp", "6" = "348-349 bp", "8" = "386-390 bp")
d[, bin_factor := factor(BIN_LABELS[as.character(bin_id)], levels = BIN_LABELS[as.character(BINS)])]
d[, sample_label := factor(SAMPLE_LABELS[sample], levels = c("CENP-A rep1", "CENP-A rep2"))]

med <- d[, .(med = median(value, na.rm = TRUE)), by = .(bin_id, bin_factor, sample, sample_label, grp)]
med_wide <- dcast(med, bin_id + bin_factor + sample + sample_label ~ grp, value.var = "med")

perm_pvals <- fread(file.path(RESULTS_DIR, "period_enrichment_merged_permutation_results.csv"))
dir <- merge(med_wide,
             perm_pvals[bin_id %in% BINS & sample %in% SAMPLES,
                        .(bin_id, sample, direction = ifelse(paired_effect > 0, "enriched", "depleted"))],
             by = c("bin_id", "sample"))
n_sig <- perm_pvals[bin_id %in% BINS & sample %in% SAMPLES & p_adj < 0.001, .N]

message("Observed vs null medians (log2 normalized, MERGED arrays; descriptive) and direction (from merged permutation test):")
print(copy(dir)[, .(bin_id, sample_label, direction, foreground = round(foreground, 3), null = round(null, 3))])
message(sprintf("Panels with P_adj < 0.001: %d / %d", n_sig, length(BINS) * length(SAMPLES)))

# ── Write the cache the plot script reads ────────────────────────────────────
# Factor columns (bin_factor/sample_label) are NOT saved — the plot script
# rebuilds them from the static label maps, so the cache holds only data.
cache <- d[, .(bin_id, sample, grp, value)]
fwrite(cache, DATA_FILE, sep = "\t")
message("Saved cache: ", DATA_FILE, " (", nrow(cache), " rows x ",
        ncol(cache), " cols)")
message("Done.")
