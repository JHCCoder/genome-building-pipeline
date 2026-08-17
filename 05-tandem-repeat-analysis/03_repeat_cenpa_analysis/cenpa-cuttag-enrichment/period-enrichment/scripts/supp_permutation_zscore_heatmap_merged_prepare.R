#!/usr/bin/env Rscript
# ============================================================================
# supp_permutation_zscore_heatmap_merged_prepare.R
# DATA-PROCESSING half of the per-chromosome x per-bin CENP-A enrichment
# Z-score heatmap (MERGED arrays, rep1 = XG_150).
#
# This script contains ALL of the expensive computation that used to live
# inside supp_permutation_zscore_heatmap_merged.R:
#
#   1. Global pseudocount from all 4 merged samples (foreground signal files).
#   2. Per-(chrom,bin) observed median (rep1) and null mean/SD from the
#      chromosome- and length-matched shuffled MERGED ARRAYS (1,000 iters):
#        - bins 3-5  -> ctrl (5,000-array subsample) null
#        - bins 6-9  -> full-set null
#   3. Per-(chrom,bin) Z-score, capped at +/-5.
#   4. Per-chromosome permutation test (star basis): matched design for ctrl
#      bins, full-set for bins 6-9; BH-FDR within each bin across chromosomes.
#
# It writes the result grid to
#     ../results/supp_permutation_zscore_heatmap_merged_matrix.tsv
#
# The PLOTTING half is scripts/supp_permutation_zscore_heatmap_merged.R, which
# reads this matrix and renders the figure in seconds. Run THIS script only
# when the raw data (signal/permutation counts, library sizes, binning)
# changes; for any pure plotting tweak just re-run the plot script.
#
# Per-chromosome test power caveat: a chromosome with only a handful of arrays
# in a sparse bin (e.g. bin 6) is tested against a null built from those same
# few shuffled arrays, so it rarely reaches significance — an honest reflection
# of low per-chromosome replication in sparse bins.
#
# Z = (observed per-chromosome median - null mean) / (null SD + 1e-10), null =
# per-iteration median of chromosome- and length-matched shuffled MERGED ARRAYS
# (1,000 iters). Bins 3-5 use the 5,000-array subsample (ctrl) null, bins 6-9
# the full-set null — same sources as the analysis. Z capped at +/-5.
#
# Why some bin-6 (348-349 bp) cells are grey: that bin has only 799 merged
# arrays genome-wide, and chr13/14/19/23 contain ZERO 348-349 bp arrays, so
# there is no observed signal and nothing to shuffle there. chrY has no TRF
# arrays in any bin 3-9 and is therefore absent from the whole heatmap.
#
# Usage (run rarely — only when raw data changes):
#   conda activate r-visualizations
#   Rscript supp_permutation_zscore_heatmap_merged_prepare.R
#
# Saved: ../results/supp_permutation_zscore_heatmap_merged_matrix.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

PERIOD_DIR    <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment"
PARENT_DIR    <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
MERGED_DIR    <- file.path(PERIOD_DIR, "data", "merged")
COUNTS_DIR    <- file.path(MERGED_DIR, "counts")
BG_COUNTS_DIR <- file.path(MERGED_DIR, "permutation", "counts")
RESULTS_DIR   <- file.path(PERIOD_DIR, "results")

S <- "XG_150"                       # heatmap is rep1 (as in the original figure)
ALL_SAMPLES  <- c("XG_150", "XG_151", "XG_152", "XG_153")
PLOT_BINS    <- 3:9                 # bins 1-2 (micro/minisatellites) removed
CTRL_BINS    <- 3:5                 # 5,000-array subsample null
FULL_BINS    <- 6:9                 # full-set null

MATRIX_FILE <- file.path(RESULTS_DIR,
                         "supp_permutation_zscore_heatmap_merged_matrix.tsv")

# ── Library sizes ─────────────────────────────────────────────────────────────
lib_file <- file.path(PARENT_DIR, "data", "qc", "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
} else {
  stop("Library size file not found: ", lib_file)
}

# ── Global pseudocount (all 4 merged samples — matches the analysis) ──────────
sig_all <- rbindlist(lapply(ALL_SAMPLES, function(s) {
  x <- fread(file.path(COUNTS_DIR, paste0("trf_signal_array_", s, ".tsv")), header = FALSE,
             col.names = c("chrom", "start", "end", "array_id", "bin_id",
                           "n_intervals", "mean_coverage"))
  x[, norm_signal := mean_coverage / (lib_sizes[s] / 1e6)]
}))
pseudo <- min(sig_all$norm_signal[sig_all$norm_signal > 0], na.rm = TRUE) / 2
message(sprintf("Global pseudocount (merged arrays): %.6e (log2 = %.3f)", pseudo, log2(pseudo)))

# ── Foreground (rep1) ─────────────────────────────────────────────────────────
fg <- fread(file.path(COUNTS_DIR, paste0("trf_signal_array_", S, ".tsv")), header = FALSE,
            col.names = c("chrom", "start", "end", "array_id", "bin_id",
                          "n_intervals", "mean_coverage"))
fg[, norm_signal := mean_coverage / (lib_sizes[S] / 1e6)]
fg[, log2_signal := log2(norm_signal + pseudo)]
chr_fg <- fg[bin_id %in% PLOT_BINS,
             .(fg_median = median(log2_signal, na.rm = TRUE)),
             by = .(chrom, bin_id)]

# ── Background (stream-filtered; never fully loaded) ─────────────────────────
# For the per-chromosome test we need the observed signal of the SUBSAMPLED
# arrays too (matched design, same as the bin-level analysis), so ctrl rows keep
# array_id and are joined to the observed log2 signal by array_id.
fg_obs <- fg[, .(array_id, obs_log2 = log2_signal)]
bg_rows <- list()
for (b in CTRL_BINS) {                       # ctrl null, WITH header
  f <- file.path(COUNTS_DIR, paste0("trf_ctrl_bg_signal_array_", S, ".tsv"))
  x <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || $1==%d' '%s'", b, f), header = TRUE)
  setnames(x, "mean_signal", "mean_coverage")
  x[, norm_signal := mean_coverage / (lib_sizes[S] / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
  x <- merge(x, fg_obs, by = "array_id", all.x = TRUE)
  bg_rows[[length(bg_rows) + 1]] <- x[, .(chrom, bin_id, iter, log2_signal, obs_log2)]
}
for (b in FULL_BINS) {                       # full-set null, NO header
  f <- file.path(BG_COUNTS_DIR, paste0("trf_bg_signal_array_", S, ".tsv"))
  x <- fread(cmd = sprintf("awk -F'\\t' '$5==%d' '%s'", b, f), header = FALSE,
             col.names = c("chrom", "start", "end", "array_id", "bin_id", "iter", "mean_coverage"))
  x[, norm_signal := mean_coverage / (lib_sizes[S] / 1e6)]
  x[, log2_signal := log2(norm_signal + pseudo)]
  bg_rows[[length(bg_rows) + 1]] <- x[, .(chrom, bin_id, iter, log2_signal)]
}
bg <- rbindlist(bg_rows, fill = TRUE)

# ── Per-(chrom,bin) Z-scores (formula identical to the analysis) ─────────────
chr_bg <- bg[, .(bg_median = median(log2_signal, na.rm = TRUE)),
             by = .(chrom, bin_id, iter)]
chr_bg_null <- chr_bg[, .(null_mean = mean(bg_median), null_sd = sd(bg_median)),
                      by = .(chrom, bin_id)]
chr_z <- merge(chr_fg, chr_bg_null, by = c("chrom", "bin_id"))
chr_z[, z_score := (fg_median - null_mean) / (null_sd + 1e-10)]
chr_z[null_sd == 0, z_score := NA_real_]
chr_z[, z_capped := ifelse(is.na(z_score), NA_real_, pmax(pmin(z_score, 5), -5))]

# ── Complete grid: all (chrom x bin) so missing cells become explicit NA ──────
# Keep the same chromosome set as the original figure (bin-7 anchor = every
# chromosome with centromeric-satellite data). chrY has no bin 3-9 arrays at
# all, so it is not present.
chr_set <- chr_z[bin_id == 7, unique(chrom)]
grid <- CJ(chrom = chr_set, bin_id = PLOT_BINS)
chr_z <- merge(grid, chr_z, by = c("chrom", "bin_id"), all.x = TRUE)
chr_z[is.na(z_capped) & is.na(fg_median), z_why := "no arrays"]   # descriptive

# ── Per-chromosome permutation test (the star basis for each cell) ────────────
# Same empirical test as the bin-level analysis, restricted to each chromosome:
#   * ctrl bins (3-5): matched — per iteration, over the SUBSAMPLED arrays of
#     chromosome c, D = median(observed) - median(shuffled) on the same arrays.
#   * full bins (6-9): shuffled set = all arrays; observed = chromosome's
#     full-bin median (constant over iterations).
#   p_enrich = mean(D <= 0); p_deplete = mean(D >= 0); two-sided p = 2*min;
#   BH-corrected WITHIN each bin across chromosomes (FDR per bin).
obs_full <- fg[bin_id %in% FULL_BINS, .(obs_full = median(log2_signal, na.rm = TRUE)),
               by = .(chrom, bin_id)]
D_list <- list()
ctrl_med <- bg[!is.na(obs_log2),
               .(obs_med = median(obs_log2, na.rm = TRUE),
                 bg_med = median(log2_signal, na.rm = TRUE)),
               by = .(chrom, bin_id, iter)]
ctrl_med <- ctrl_med[!is.na(obs_med) & !is.na(bg_med)]
D_list[["ctrl"]] <- ctrl_med[, .(chrom, bin_id, iter, D = obs_med - bg_med)]
full_med <- bg[is.na(obs_log2),
               .(bg_med = median(log2_signal, na.rm = TRUE)),
               by = .(chrom, bin_id, iter)]
full_med <- merge(full_med, obs_full, by = c("chrom", "bin_id"))
full_med <- full_med[!is.na(bg_med)]
D_list[["full"]] <- full_med[, .(chrom, bin_id, iter, D = obs_full - bg_med)]
allD <- rbindlist(D_list)

pchr <- allD[, .(p_enrich = mean(D <= 0), p_deplete = mean(D >= 0), n_iters = .N),
             by = .(chrom, bin_id)]
pchr[, p_two := pmin(p_enrich, p_deplete) * 2]
pchr[, p_two := pmin(p_two, 1)]
pchr[, p_adj := p.adjust(p_two, method = "BH"), by = bin_id]
pchr[, sig := ifelse(p_adj < 0.001, "***", ifelse(p_adj < 0.01, "**",
                     ifelse(p_adj < 0.05, "*", "ns")))]

message("Per-chromosome significance (rep1; BH FDR within each bin):")
sig_tab <- pchr[sig != "ns", .(n_sig = .N), by = bin_id][order(bin_id)]
print(sig_tab, row.names = FALSE)
message("Effective iterations per (chrom,bin) cell: range",
        paste(range(allD[, .(n_iters = .N), by = .(chrom, bin_id)]$n_iters), collapse = " - "))

# Attach per-chromosome significance to the plotting grid
chr_z <- merge(chr_z, pchr[, .(chrom, bin_id, p_adj, sig, n_iters)],
               by = c("chrom", "bin_id"), all.x = TRUE)
chr_z[, sig_mark := ifelse(is.na(sig) | sig == "ns", "", sig)]

# ── Write the cache the plot script reads ─────────────────────────────────────
# Factors (chrom_factor/bin_factor) are NOT saved — the plot script rebuilds
# them from chrom/bin_id, so the cache holds only data, not view state.
cache <- chr_z[, .(chrom, bin_id, fg_median, null_mean, null_sd,
                   z_score, z_capped, z_why, p_adj, sig, n_iters, sig_mark)]
fwrite(cache, MATRIX_FILE, sep = "\t")
message("Saved matrix: ", MATRIX_FILE, " (", nrow(cache), " rows x ",
        ncol(cache), " cols)")
message("Done.")
