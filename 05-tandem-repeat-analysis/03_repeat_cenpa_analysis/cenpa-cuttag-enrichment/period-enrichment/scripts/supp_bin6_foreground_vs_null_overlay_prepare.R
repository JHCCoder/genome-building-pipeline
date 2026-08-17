#!/usr/bin/env Rscript
# ============================================================================
# supp_bin6_foreground_vs_null_overlay_prepare.R
# DATA-PROCESSING half of the per-interval overlay for the 348-349 bp monomer
# (bin 6):
#   * grey  = NULL per-interval signal  (all shuffled intervals across 1,000
#             iterations, pooled; data/permutation/per_bin/full/bin6_*.tsv)
#   * color = observed CENP-A foreground per-interval signal
#             (data/counts/trf_signal_*.tsv, bin 6 only)
#
# This script contains ALL of the expensive computation that used to live
# inside supp_bin6_foreground_vs_null_overlay.R: the global pseudocount from
# all 4 samples, the foreground and null per-interval log2 signals for bin 6
# (CENP-A samples), and the median printouts. It writes the plotting data to
#     ../results/supp_bin6_foreground_vs_null_distribution_data.tsv
#     columns: grp ("foreground"/"null"), sample, log2_signal
#
# The PLOTTING half is scripts/supp_bin6_foreground_vs_null_overlay.R, which
# reads this table (plus the small permutation results CSV for the P-value
# annotation) and renders the figure in seconds. Run THIS script only when the
# raw data changes; for any pure plotting tweak just re-run the plot script.
#
# Normalization matches 3_analyze_period_enrichment.R EXACTLY (verified to
# machine precision against results/period_enrichment_permutation_results.csv):
#   norm_signal = mean_coverage / (lib_size / 1e6)
#   pseudo      = min(norm_signal > 0) / 2  (global, all 4 samples)
#   log2_signal = log2(norm_signal + pseudo)
#
# Usage (run rarely — only when raw data changes):
#   conda activate r-visualizations
#   Rscript supp_bin6_foreground_vs_null_overlay_prepare.R
#
# Saved: ../results/supp_bin6_foreground_vs_null_distribution_data.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

PERIOD_DIR  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PARENT_DIR  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
COUNTS_DIR  <- file.path(PERIOD_DIR, "data", "counts")
PER_BIN_DIR <- file.path(PERIOD_DIR, "data", "permutation", "per_bin", "full")
RESULTS_DIR <- file.path(PERIOD_DIR, "results")

SAMPLES <- c("XG_150", "XG_151")                       # CENP-A replicates
ALL_SAMPLES   <- c("XG_150", "XG_151", "XG_152", "XG_153")
BIN6 <- 6  # 348-349 bp = the 349-bp monomer family in the 9-bin scheme

DATA_FILE <- file.path(RESULTS_DIR,
                       "supp_bin6_foreground_vs_null_distribution_data.tsv")

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
sig_all[, lib_size := lib_sizes[sample]]
sig_all[, norm_signal := mean_coverage / (lib_size / 1e6)]
pseudo <- min(sig_all$norm_signal[sig_all$norm_signal > 0], na.rm = TRUE) / 2
message(sprintf("Global pseudocount: %.6e (log2 = %.3f)", pseudo, log2(pseudo)))

# ── Foreground: bin-6 per-interval signal, CENP-A samples only ─────────────
fg_bin6 <- sig_all[bin_id == BIN6 & sample %in% SAMPLES]
fg_bin6[, log2_signal := log2(norm_signal + pseudo)]

# ── Null: bin-6 per-interval signal pooled across all 1,000 iterations ─────
nul_bin6 <- rbindlist(lapply(SAMPLES, function(s) {
  f <- file.path(PER_BIN_DIR, sprintf("bin%d_%s.tsv", BIN6, s))
  x <- fread(f, header = FALSE,
             col.names = c("chrom", "start", "end", "bin_id", "iter", "mean_coverage"))
  x[, sample := s]
  x[, lib_size := lib_sizes[sample]]
  x[, norm_signal := mean_coverage / (lib_size / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
}))

message(sprintf("Foreground bin-6 intervals/sample: %d", nrow(fg_bin6) / length(SAMPLES)))
message(sprintf("Null bin-6 rows (interval x iteration) per sample: %d",
                nrow(nul_bin6) / length(SAMPLES)))

# ── Combined plotting data (one file the plot script reads) ─────────────────
d <- rbind(
  fg_bin6[, .(grp = "foreground", sample, log2_signal)],
  nul_bin6[, .(grp = "null", sample, log2_signal)]
)

# ── Medians for the vertical markers (same printout as the combined script) ──
SAMPLE_LABELS <- c("XG_150" = "CENP-A rep1", "XG_151" = "CENP-A rep2")
med <- d[, .(median = median(log2_signal, na.rm = TRUE)),
         by = .(sample, grp)]
med[, sample_label := factor(SAMPLE_LABELS[sample], levels = c("CENP-A rep1", "CENP-A rep2"))]
message("Foreground medians:")
print(med[grp == "foreground", .(sample, sample_label, median)])
message("Null medians:")
print(med[grp == "null", .(sample, sample_label, median)])

# ── Write the cache the plot script reads ────────────────────────────────────
# Factor columns are NOT saved — the plot script rebuilds them from the static
# label maps, so the cache holds only data.
fwrite(d, DATA_FILE, sep = "\t")
message("Saved cache: ", DATA_FILE, " (", nrow(d), " rows x ",
        ncol(d), " cols)")
message("Done.")
