#!/usr/bin/env Rscript
# ============================================================================
# supp_monomer_bins_foreground_vs_null_overlay_prepare.R
# DATA-PROCESSING half of the per-interval overlay for the three CENP-A monomer
# families: bins 4 (193-195 bp), 6 (348-349 bp) and 8 (386-390 bp).
#
# This script contains ALL of the expensive computation that used to live
# inside supp_monomer_bins_foreground_vs_null_overlay.R:
#   * global pseudocount from all 4 samples (foreground signal files)
#   * foreground per-interval log2 signal for bins 4/6/8, CENP-A samples
#   * null per-interval log2 signal (pooled across 1,000 shuffled iterations)
#     from the ctrl/full per-bin permutation files; nulls >1M rows thinned to
#     1M (set.seed 42) — the thinning is part of the cached data.
#   * descriptive observed vs null medians and panel direction/n_sig printouts
#     (direction read from the PERMUTATION results CSV).
#
# It writes the plotting data to
#     ../results/supp_monomer_bins_foreground_vs_null_distribution_data.tsv
#     columns: bin_id, sample, grp ("foreground"/"null"), value (log2 normalized)
#
# The PLOTTING half is scripts/supp_monomer_bins_foreground_vs_null_overlay.R,
# which reads this table (plus the small permutation results CSV for the panel
# labels) and renders the figure in seconds. Run THIS script only when the raw
# data changes; for any pure plotting tweak just re-run the plot script.
#
# Normalization matches 3_analyze_period_enrichment.R EXACTLY (verified to
# machine precision against results/period_enrichment_permutation_results.csv):
#   norm_signal = mean_coverage / (lib_size / 1e6)
#   pseudo      = min(norm_signal > 0) / 2  (global, all 4 samples)
#   log2_signal = log2(norm_signal + pseudo)
#
# NOTE on the zero-coverage spike at log2 = -9.38: a real property of the
# shuffled background (matches bg_mean_of_medians in the permutation CSV), not
# an artifact of this figure.
#
# Usage (run rarely — only when raw data changes):
#   conda activate r-visualizations
#   Rscript supp_monomer_bins_foreground_vs_null_overlay_prepare.R
#
# Saved: ../results/supp_monomer_bins_foreground_vs_null_distribution_data.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

PERIOD_DIR  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PARENT_DIR  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
COUNTS_DIR  <- file.path(PERIOD_DIR, "data", "counts")
PER_BIN_DIR <- file.path(PERIOD_DIR, "data", "permutation", "per_bin")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

SAMPLES <- c("XG_150", "XG_151")                       # CENP-A replicates
ALL_SAMPLES   <- c("XG_150", "XG_151", "XG_152", "XG_153")
BINS <- c(4, 6, 8)   # 193-195 bp, 348-349 bp, 386-390 bp

DATA_FILE <- file.path(RESULTS_DIR,
                       "supp_monomer_bins_foreground_vs_null_distribution_data.tsv")

# ── Library sizes (must match the analysis) ────────────────────────────────
lib_file <- file.path(PARENT_DIR, "data", "qc", "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
} else {
  stop("Library size file not found: ", lib_file)
}

# ── Global pseudocount (all samples, all bins — same as the analysis) ──────
sig_all <- rbindlist(lapply(ALL_SAMPLES, function(s) {
  x <- fread(file.path(COUNTS_DIR, paste0("trf_signal_", s, ".tsv")), header = FALSE,
             col.names = c("chrom", "start", "end", "interval_id", "bin_id",
                           "period_size", "copies_aligned", "match_percent",
                           "mean_coverage"))
  x[, sample := s]
}))
sig_all[, norm_signal := mean_coverage / (lib_sizes[sample] / 1e6)]
pseudo <- min(sig_all$norm_signal[sig_all$norm_signal > 0], na.rm = TRUE) / 2
message(sprintf("Global pseudocount: %.6e (log2 = %.3f)", pseudo, log2(pseudo)))

lgsig <- function(cov, s) log2(cov / (lib_sizes[s] / 1e6) + pseudo)

# ── Foreground: per-interval signal for bins 4/6/8, CENP-A samples ──────────
fg_all <- sig_all[bin_id %in% BINS & sample %in% SAMPLES]
fg_all[, value := lgsig(mean_coverage, sample)]

# ── Null per-interval (pooled across 1,000 iterations) ─────────────────────
load_null <- function(s, b) {
  if (b %in% 1:5) {  # bins 1-5: ctrl (subsample) null, 7 cols WITH header
    x <- fread(file.path(PER_BIN_DIR, "ctrl", sprintf("bin%d_%s.tsv", b, s)), header = TRUE)
    setnames(x, "mean_signal", "mean_coverage")
  } else {           # bins 6-9: full-set null, 6 cols NO header
    x <- fread(file.path(PER_BIN_DIR, "full", sprintf("bin%d_%s.tsv", b, s)), header = FALSE,
               col.names = c("chrom", "start", "end", "bin_id", "iter", "mean_coverage"))
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
BIN_LABELS <- c("4" = "193-195 bp (bin 4)", "6" = "348-349 bp (bin 6)", "8" = "386-390 bp (bin 8)")
d[, bin_factor := factor(BIN_LABELS[as.character(bin_id)], levels = BIN_LABELS[as.character(BINS)])]
d[, sample_label := factor(SAMPLE_LABELS[sample], levels = c("CENP-A rep1", "CENP-A rep2"))]

med <- d[, .(med = median(value, na.rm = TRUE)), by = .(bin_id, bin_factor, sample, sample_label, grp)]
med_wide <- dcast(med, bin_id + bin_factor + sample + sample_label ~ grp, value.var = "med")

perm_pvals <- fread(file.path(RESULTS_DIR, "period_enrichment_permutation_results.csv"))
dir <- merge(med_wide,
             perm_pvals[bin_id %in% BINS & sample %in% SAMPLES,
                        .(bin_id, sample, direction = ifelse(paired_effect > 0, "enriched", "depleted"))],
             by = c("bin_id", "sample"))
n_sig <- perm_pvals[bin_id %in% BINS & sample %in% SAMPLES & p_adj < 0.001, .N]

message("Observed vs null medians (log2 normalized; descriptive) and direction (from permutation test):")
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
