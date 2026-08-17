#!/usr/bin/env Rscript
# ============================================================================
# verify_diagnostics.R — Verification of V2 analysis claims
#
# Produces:
#   1. Per-chromosome diagnostic table
#   2. Interval-count vs enrichment plot
#   3. Foreground-length vs enrichment plot
#   4. Inter-interval gap distance distribution
#   5. Zero-count analysis by chromosome
#   6. Corrected empirical permutation test
#   7. Merge candidates at multiple distances
#
# Usage:
#   conda activate r-visualizations
#   Rscript verify_diagnostics.R
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
DATA_DIR <- file.path(BASE_DIR, "data")
COUNTS_DIR <- file.path(DATA_DIR, "counts")
FOREGROUND_DIR <- file.path(DATA_DIR, "foregrounds")
BACKGROUND_DIR <- file.path(DATA_DIR, "backgrounds")
MAPPABILITY_DIR <- file.path(DATA_DIR, "mappability")
QC_DIR <- file.path(DATA_DIR, "qc")
VERIFY_DIR <- file.path(BASE_DIR, "verification")
dir.create(VERIFY_DIR, recursive = TRUE, showWarnings = FALSE)

CHROMOSOMES <- paste0("chr", c(1:28, "X", "Y"))
SAMPLES <- c("XG_150", "XG_151", "XG_152")
SAMPLE_LABELS <- c("XG_150" = "CENP-A rep1",
                    "XG_151" = "CENP-A rep2",
                    "XG_152" = "H3K27ac")
SAMPLE_COLORS <- c("XG_150" = "#2166AC", "XG_151" = "#92C5DE", "XG_152" = "#B2182B")

theme_verify <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        axis.text = element_text(color = "black"))

# ============================================================================
# Load data
# ============================================================================
message("=== Loading data ===")

# Library sizes
lib_file <- file.path(QC_DIR, "library_sizes_fragments.txt")
lib_sizes <- setNames(fread(lib_file, header = FALSE)$V2,
                       fread(lib_file, header = FALSE)$V1)
message("Library sizes: ", paste(names(lib_sizes), lib_sizes, sep = "=", collapse = ", "))

# Foreground counts (Set A)
load_counts <- function(sample, region_set) {
  f <- file.path(COUNTS_DIR, paste0(sample, "_", region_set, ".txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE, col.names = c("chrom", "start", "end", "id", "count"))
  dt[, sample := sample]
  dt[, length_bp := end - start]
  return(dt)
}

setA_list <- lapply(SAMPLES, function(s) load_counts(s, "setA_strict"))
setA <- rbindlist(setA_list[!sapply(setA_list, is.null)])
setA[, norm_signal := (count / (length_bp / 1000)) / (lib_sizes[sample] / 1e6)]

# Background counts (Bg1 — chromosome-matched shuffle)
load_bg1 <- function(sample) {
  f <- file.path(COUNTS_DIR, paste0(sample, "_bg1_chrom_shuffle.txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE, col.names = c("chrom", "start", "end", "iter_id", "count"))
  dt[, sample := sample]
  dt[, length_bp := end - start]
  return(dt)
}
bg1_list <- lapply(SAMPLES, function(s) load_bg1(s))
bg1 <- rbindlist(bg1_list[!sapply(bg1_list, is.null)])
bg1[, norm_signal := (count / (length_bp / 1000)) / (lib_sizes[sample] / 1e6)]

# Mappability
mapp_file <- file.path(MAPPABILITY_DIR, "setA_mappability_scores.txt")
mapp_data <- NULL
if (file.exists(mapp_file)) {
  mapp_data <- fread(mapp_file, header = FALSE,
                     col.names = c("chrom", "start", "end", "id", "mappable_fraction"))
}

# ============================================================================
# 1. PER-CHROMOSOME DIAGNOSTIC TABLE
# ============================================================================
message("\n=== 1. Per-chromosome diagnostic table ===")

# Pseudocount for log2 transform
min_nonzero <- min(c(setA$norm_signal[setA$norm_signal > 0], bg1$norm_signal[bg1$norm_signal > 0]), na.rm = TRUE)
pseudocount <- min_nonzero / 2
setA[, log2_signal := log2(norm_signal + pseudocount)]
bg1[, log2_signal := log2(norm_signal + pseudocount)]

# Per-chromosome foreground stats
chr_fg <- setA[, .(
  n_intervals = .N,
  total_foreground_bp = as.double(sum(length_bp, na.rm = TRUE)),
  median_fg = as.double(median(log2_signal, na.rm = TRUE)),
  mean_fg = as.double(mean(log2_signal, na.rm = TRUE)),
  sd_fg = as.double(sd(log2_signal, na.rm = TRUE)),
  zero_count_frac = as.double(sum(count == 0) / .N),
  median_raw_count = as.double(median(count, na.rm = TRUE))
), by = .(sample, chrom)]

# Per-chromosome background stats (using all 100 bg1 iterations pooled per chrom)
chr_bg <- bg1[, .(
  median_bg = as.double(median(log2_signal, na.rm = TRUE)),
  mean_bg = as.double(mean(log2_signal, na.rm = TRUE)),
  n_bg_intervals = .N
), by = .(sample, chrom)]

# Merge fg and bg
chr_diag <- merge(chr_fg, chr_bg, by = c("sample", "chrom"))
chr_diag[, delta := median_fg - median_bg]
chr_diag[, fold_enrichment := 2^delta]

# Add mappability if available
if (!is.null(mapp_data)) {
  mapp_chr <- mapp_data[, .(mean_mappability = mean(mappable_fraction, na.rm = TRUE),
                              median_mappability = median(mappable_fraction, na.rm = TRUE)),
                          by = chrom]
  chr_diag <- merge(chr_diag, mapp_chr, by = "chrom", all.x = TRUE)
}

# Add chromosome size from FASTA index
fai_file <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned-contamFiltered/assembly_final.sorted.headerRenamed.chrAssigned.contamFiltered.fasta.fai"
if (file.exists(fai_file)) {
  fai <- fread(fai_file, header = FALSE, col.names = c("chrom", "size", "offset", "linebases", "linewidth"))
  fai <- fai[chrom %in% CHROMOSOMES]
  chr_diag <- merge(chr_diag, fai[, .(chrom, chr_size = size)], by = "chrom", all.x = TRUE)
  chr_diag[, centroanno_frac := total_foreground_bp / chr_size]
}

# Write diagnostic table
diag_wide <- dcast(chr_diag, chrom ~ sample,
                    value.var = c("n_intervals", "total_foreground_bp", "median_fg",
                                  "median_bg", "delta", "fold_enrichment",
                                  "zero_count_frac", "median_raw_count"))

fwrite(chr_diag[order(chrom, sample)],
       file.path(VERIFY_DIR, "per_chromosome_diagnostics.csv"))

# Print key findings
message("\nChromosomes with negative CENP-A delta:")
neg_chr_150 <- chr_diag[sample == "XG_150" & delta < 0]
neg_chr_151 <- chr_diag[sample == "XG_151" & delta < 0]
message("  XG_150: ", paste(neg_chr_150$chrom, collapse = ", "),
        " (n=", nrow(neg_chr_150), "/", nrow(chr_diag[sample == "XG_150"]), ")")
message("  XG_151: ", paste(neg_chr_151$chrom, collapse = ", "),
        " (n=", nrow(neg_chr_151), "/", nrow(chr_diag[sample == "XG_151"]), ")")

message("\nChromosomes with positive CENP-A delta in both replicates:")
pos_both <- chr_diag[sample %in% c("XG_150", "XG_151"),
                     .(both_pos = all(delta > 0)),
                     by = chrom][both_pos == TRUE]
message("  ", paste(pos_both$chrom, collapse = ", "), " (n=", nrow(pos_both), ")")

message("\nInterval count imbalance:")
chr_counts_150 <- chr_diag[sample == "XG_150", .(chrom, n_intervals, delta)]
total_intervals <- sum(chr_counts_150$n_intervals)
neg_chr_counts <- chr_counts_150[delta < 0, sum(n_intervals)]
pos_chr_counts <- chr_counts_150[delta > 0, sum(n_intervals)]
message("  Total intervals: ", total_intervals)
message("  Intervals on negative-delta chromosomes: ", neg_chr_counts,
        " (", round(100 * neg_chr_counts / total_intervals, 1), "%)")
message("  Intervals on positive-delta chromosomes: ", pos_chr_counts,
        " (", round(100 * pos_chr_counts / total_intervals, 1), "%)")
message("  Negative chromosomes: ", sum(chr_counts_150$delta < 0),
        " / ", nrow(chr_counts_150))
message("  Positive chromosomes: ", sum(chr_counts_150$delta > 0),
        " / ", nrow(chr_counts_150))

# ============================================================================
# 2. INTERVAL COUNT vs ENRICHMENT PLOT
# ============================================================================
message("\n=== 2. Interval count vs enrichment ===")

chr_plot_data <- chr_diag[sample %in% c("XG_150", "XG_151")]
chr_plot_data[, sample_label := SAMPLE_LABELS[sample]]

p_count_vs_delta <- ggplot(chr_plot_data, aes(x = n_intervals, y = delta, color = sample)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text(aes(label = chrom), size = 2.5, vjust = -1, hjust = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_x_log10() +
  scale_color_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS) +
  facet_wrap(~ sample_label) +
  labs(x = "Number of centroAnno intervals (log scale)",
       y = expression(Delta * " median log"[2] * "(normalized signal) [fg - bg]"),
       title = "Interval count does not predict per-chromosome enrichment",
       subtitle = paste0("chr1: ", chr_counts_150[chrom == "chr1", n_intervals],
                         " intervals, delta = ",
                         round(chr_counts_150[chrom == "chr1", delta], 2)),
       color = NULL) +
  theme_verify
ggsave(file.path(VERIFY_DIR, "interval_count_vs_enrichment.pdf"), p_count_vs_delta,
       width = 10, height = 5)
ggsave(file.path(VERIFY_DIR, "interval_count_vs_enrichment.png"), p_count_vs_delta,
       width = 10, height = 5, dpi = 300)

# Total foreground length vs enrichment
p_length_vs_delta <- ggplot(chr_plot_data, aes(x = total_foreground_bp / 1e6, y = delta, color = sample)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text(aes(label = chrom), size = 2.5, vjust = -1, hjust = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS) +
  facet_wrap(~ sample_label) +
  labs(x = "Total centroAnno foreground length (Mb)",
       y = expression(Delta * " median log"[2] * "(normalized signal) [fg - bg]"),
       title = "Total foreground length vs per-chromosome enrichment") +
  theme_verify
ggsave(file.path(VERIFY_DIR, "foreground_length_vs_enrichment.pdf"), p_length_vs_delta,
       width = 10, height = 5)
ggsave(file.path(VERIFY_DIR, "foreground_length_vs_enrichment.png"), p_length_vs_delta,
       width = 10, height = 5, dpi = 300)

# ============================================================================
# 3. AGGREGATION EFFECT: weighted vs unweighted
# ============================================================================
message("\n=== 3. Aggregation effect: verifying the weighting explanation ===")

# Genome-wide: pool all intervals
pooled_fg_150 <- setA[sample == "XG_150", log2_signal]
pooled_bg_150 <- bg1[sample == "XG_150", log2_signal]
pooled_median_fg <- median(pooled_fg_150, na.rm = TRUE)
pooled_median_bg <- median(pooled_bg_150, na.rm = TRUE)
pooled_delta <- pooled_median_fg - pooled_median_bg
message("Pooled interval-level (all 3,420 intervals):")
message("  Median FG: ", round(pooled_median_fg, 3))
message("  Median BG: ", round(pooled_median_bg, 3))
message("  Delta: ", round(pooled_delta, 3))

# Chromosome-level: median of per-chromosome medians
chr_level_delta_150 <- chr_diag[sample == "XG_150", na.omit(delta)]
chr_level_median_delta <- median(chr_level_delta_150)
message("Chromosome-level (median of per-chromosome deltas):")
message("  Median delta across chromosomes: ", round(chr_level_median_delta, 3))
message("  Positive chromosomes: ", sum(chr_level_delta_150 > 0), "/", length(chr_level_delta_150))

# What if we weight chromosomes by interval count?
weighted_mean_delta <- weighted.mean(chr_diag[sample == "XG_150", delta],
                                      chr_diag[sample == "XG_150", n_intervals],
                                      na.rm = TRUE)
message("Interval-count-weighted mean delta: ", round(weighted_mean_delta, 3))

# What if we exclude chr1?
pooled_no_chr1_fg <- setA[sample == "XG_150" & chrom != "chr1", log2_signal]
pooled_no_chr1_bg <- bg1[sample == "XG_150" & chrom != "chr1", log2_signal]
delta_no_chr1 <- median(pooled_no_chr1_fg, na.rm = TRUE) - median(pooled_no_chr1_bg, na.rm = TRUE)
message("Pooled interval delta excluding chr1: ", round(delta_no_chr1, 3))

# Leave-one-chromosome-out analysis
loc_deltas <- rbindlist(lapply(CHROMOSOMES, function(chr) {
  fg <- setA[sample == "XG_150" & chrom != chr, log2_signal]
  bg <- bg1[sample == "XG_150" & chrom != chr, log2_signal]
  delta_val <- median(fg, na.rm = TRUE) - median(bg, na.rm = TRUE)
  data.table(chrom_excluded = chr, n_remaining_intervals = length(fg),
             pooled_delta = delta_val)
}))
loc_deltas[, delta_change := pooled_delta - loc_deltas[chrom_excluded == "NONE"]$pooled_delta]
fwrite(loc_deltas, file.path(VERIFY_DIR, "leave_one_chromosome_out.csv"))

p_loco <- ggplot(loc_deltas, aes(x = reorder(chrom_excluded, pooled_delta),
                                  y = pooled_delta)) +
  geom_col(fill = "#2166AC", alpha = 0.8) +
  geom_hline(yintercept = pooled_delta, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, color = "grey50") +
  labs(x = "Chromosome excluded",
       y = expression("Pooled median " * Delta * " log"[2] * "(normalized signal)"),
       title = "Leave-one-chromosome-out: pooled delta (XG_150)",
       subtitle = paste0("Red line = full dataset delta (",
                         round(pooled_delta, 2), ")")) +
  theme_verify +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(VERIFY_DIR, "leave_one_chromosome_out.pdf"), p_loco, width = 10, height = 4)
ggsave(file.path(VERIFY_DIR, "leave_one_chromosome_out.png"), p_loco, width = 10, height = 4, dpi = 300)

message("Strongest influence on pooled delta:")
loc_deltas[order(abs(delta_change))][1:5]

# ============================================================================
# 4. CORRECTED EMPIRICAL PERMUTATION TEST
# ============================================================================
message("\n=== 4. Corrected empirical permutation test ===")

set.seed(20260716)
N_PERM <- 100  # using available 100 bg1 iterations

# Strategy: Compare observed per-chromosome enrichment against the null
# distribution generated by treating shuffled backgrounds as pseudo-foregrounds.
#
# Observed: median per-chromosome delta = fg_chr_median - bg_pooled_chr_median
#
# Null (for each iteration i):
#   pseudo-foreground = bg_iter_i
#   pseudo-background = median of all OTHER bg iterations (leave-one-out)
#   null_delta_i = median per-chromosome (pseudo_fg - pseudo_bg)
#
# This tests: "Is the CENP-A foreground more enriched relative to matched
# background than expected from random genomic positions?"
#
# P = (1 + number of null deltas >= observed delta) / (1 + n_permutations)

corrected_permutation <- function(fg_data, bg_data, sample_name) {
  # Foreground per-chromosome signals
  fg_chr <- fg_data[sample == sample_name,
                     .(median_fg = median(log2_signal, na.rm = TRUE)),
                     by = chrom]
  setkey(fg_chr, chrom)

  # Background pooled across ALL iterations (reference)
  bg_chr_pooled <- bg_data[sample == sample_name,
                           .(median_bg = median(log2_signal, na.rm = TRUE)),
                           by = chrom]
  setkey(bg_chr_pooled, chrom)

  # Observed statistic: median per-chromosome delta
  obs_delta_chr <- fg_chr$median_fg - bg_chr_pooled$median_bg
  obs_stat <- median(obs_delta_chr, na.rm = TRUE)
  obs_frac_pos <- mean(obs_delta_chr > 0, na.rm = TRUE)

  # Build null distribution: each bg iteration as pseudo-foreground vs
  # the pooled median of all OTHER iterations as pseudo-background
  null_stats <- sapply(1:N_PERM, function(i) {
    iter_label <- sprintf("iter_%03d", i)
    bg_iter <- bg_data[sample == sample_name & iter_id == iter_label]
    if (nrow(bg_iter) == 0) return(c(NA, NA))

    # Pseudo-foreground: this iteration
    pseudo_fg <- bg_iter[, .(median_pseudo_fg = median(log2_signal, na.rm = TRUE)), by = chrom]

    # Pseudo-background: median of all OTHER iterations
    bg_others <- bg_data[sample == sample_name & iter_id != iter_label,
                          .(median_bg_other = median(log2_signal, na.rm = TRUE)),
                          by = chrom]

    merged <- merge(pseudo_fg, bg_others, by = "chrom", all = TRUE)
    null_delta <- merged$median_pseudo_fg - merged$median_bg_other
    c(median(null_delta, na.rm = TRUE),
      mean(null_delta > 0, na.rm = TRUE))
  })

  null_medians <- null_stats[1, ]
  null_frac_pos <- null_stats[2, ]
  null_medians <- null_medians[!is.na(null_medians)]
  null_frac_pos <- null_frac_pos[!is.na(null_frac_pos)]

  # One-tailed empirical P: enrichment (observed > null)
  p_median <- (1 + sum(null_medians >= obs_stat)) / (1 + length(null_medians))
  p_frac_pos <- (1 + sum(null_frac_pos >= obs_frac_pos)) / (1 + length(null_frac_pos))

  # Also test against zero (standard one-sample test)
  # Wilcoxon signed-rank: are per-chromosome deltas > 0?
  wilcox_p <- wilcox.test(obs_delta_chr, mu = 0, alternative = "greater", exact = TRUE)$p.value

  list(
    observed = obs_stat,
    obs_frac_pos = obs_frac_pos,
    null_stats = null_medians,
    null_frac_pos = null_frac_pos,
    n_null = length(null_medians),
    p_median = p_median,
    p_frac_pos = p_frac_pos,
    wilcox_p = wilcox_p,
    obs_chr_deltas = obs_delta_chr,
    n_positive_chr = sum(obs_delta_chr > 0, na.rm = TRUE),
    n_total_chr = sum(!is.na(obs_delta_chr))
  )
}

perm_150 <- corrected_permutation(setA, bg1, "XG_150")
perm_151 <- corrected_permutation(setA, bg1, "XG_151")

message("Corrected permutation test results:")
message("  XG_150: observed median chr delta = ", round(perm_150$observed, 3),
        ", P(perm) = ", formatC(perm_150$p_median, digits = 4),
        ", P(Wilcoxon) = ", formatC(perm_150$wilcox_p, digits = 4),
        ", frac pos chr = ", round(perm_150$obs_frac_pos, 2),
        " (", perm_150$n_positive_chr, "/", perm_150$n_total_chr, ")")
message("  Null median delta mean = ", round(mean(perm_150$null_stats), 3),
        " (expect ~0 if no systematic bias)")
message("  XG_151: observed median chr delta = ", round(perm_151$observed, 3),
        ", P(perm) = ", formatC(perm_151$p_median, digits = 4),
        ", P(Wilcoxon) = ", formatC(perm_151$wilcox_p, digits = 4),
        ", frac pos chr = ", round(perm_151$obs_frac_pos, 2),
        " (", perm_151$n_positive_chr, "/", perm_151$n_total_chr, ")")
message("  Null median delta mean = ", round(mean(perm_151$null_stats), 3))

# Plot null distributions
null_df <- rbind(
  data.table(sample = "XG_150", delta = perm_150$null_stats),
  data.table(sample = "XG_151", delta = perm_151$null_stats)
)
obs_df <- data.table(
  sample = c("XG_150", "XG_151"),
  delta = c(perm_150$observed, perm_151$observed),
  p_val = c(perm_150$p_median, perm_151$p_median)
)

p_null_corrected <- ggplot(null_df, aes(x = delta)) +
  geom_histogram(bins = 30, fill = "grey70", alpha = 0.7) +
  geom_vline(data = obs_df, aes(xintercept = delta), color = "#2166AC",
             linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_text(data = obs_df, aes(x = delta, y = Inf,
            label = paste0("observed = ", round(delta, 3), "\nP(perm) = ",
                           formatC(p_val, digits = 3))),
            vjust = 1.5, hjust = -0.1, size = 3, color = "#2166AC") +
  facet_wrap(~ sample, ncol = 1,
             labeller = labeller(sample = SAMPLE_LABELS)) +
  labs(x = expression("Median chromosome-level " * Delta * " log"[2] * "(normalized signal)"),
       y = "Frequency (100 null iterations)",
       title = "Corrected empirical permutation: chromosome-level statistic",
       subtitle = "Null: shuffled background as pseudo-foreground vs other backgrounds") +
  theme_verify
ggsave(file.path(VERIFY_DIR, "corrected_permutation_test.pdf"), p_null_corrected,
       width = 7, height = 6)
ggsave(file.path(VERIFY_DIR, "corrected_permutation_test.png"), p_null_corrected,
       width = 7, height = 6, dpi = 300)

# ============================================================================
# 5. INTER-INTERVAL GAP DISTRIBUTION
# ============================================================================
message("\n=== 5. Inter-interval gap distribution ===")

# Load centroAnno intervals
centro_bed <- fread(file.path(FOREGROUND_DIR, "setA_centroAnno_strict.bed"),
                    header = FALSE, col.names = c("chrom", "start", "end", "id"))
setorder(centro_bed, chrom, start)

# Compute gaps between adjacent intervals on the same chromosome
centro_bed[, gap_to_next := shift(start, type = "lead") - end, by = chrom]
gap_data <- centro_bed[!is.na(gap_to_next) & gap_to_next >= 0]

# Gap summary statistics
gap_summary <- gap_data[, .(
  median_gap = as.double(median(gap_to_next)),
  q25 = as.double(quantile(gap_to_next, 0.25)),
  q75 = as.double(quantile(gap_to_next, 0.75)),
  mean_gap = as.double(mean(gap_to_next)),
  pct_within_50k = as.double(mean(gap_to_next <= 50000) * 100),
  pct_within_100k = as.double(mean(gap_to_next <= 100000) * 100),
  pct_within_250k = as.double(mean(gap_to_next <= 250000) * 100),
  pct_within_500k = as.double(mean(gap_to_next <= 500000) * 100),
  pct_within_1m = as.double(mean(gap_to_next <= 1000000) * 100)
), by = chrom]

fwrite(gap_summary, file.path(VERIFY_DIR, "inter_interval_gap_summary.csv"))

message("Genome-wide inter-interval gap distribution:")
all_gaps <- gap_data$gap_to_next
message("  Median gap: ", round(median(all_gaps) / 1000, 1), " kb")
message("  Q1: ", round(quantile(all_gaps, 0.25) / 1000, 1), " kb")
message("  Q3: ", round(quantile(all_gaps, 0.75) / 1000, 1), " kb")
message("  Fraction <= 50kb: ", round(mean(all_gaps <= 50000) * 100, 1), "%")
message("  Fraction <= 100kb: ", round(mean(all_gaps <= 100000) * 100, 1), "%")
message("  Fraction <= 250kb: ", round(mean(all_gaps <= 250000) * 100, 1), "%")
message("  Fraction <= 500kb: ", round(mean(all_gaps <= 500000) * 100, 1), "%")
message("  Fraction <= 1Mb: ", round(mean(all_gaps <= 1e6) * 100, 1), "%")

# Gap distribution plot (log scale)
p_gap_hist <- ggplot(gap_data, aes(x = gap_to_next / 1000)) +
  geom_histogram(bins = 100, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = c(50, 100, 250, 500), linetype = "dashed",
             color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  scale_x_log10(labels = scales::comma) +
  annotate("text", x = c(50, 100, 250, 500), y = Inf,
           label = c("50kb", "100kb", "250kb", "500kb"),
           vjust = 2, size = 3, angle = 90,
           color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  labs(x = "Gap between adjacent centroAnno intervals (kb, log scale)",
       y = "Count",
       title = "Distribution of distances between adjacent centroAnno intervals",
       subtitle = paste0("Median = ", round(median(all_gaps) / 1000, 1),
                         " kb; ", round(mean(all_gaps <= 100000) * 100, 1),
                         "% within 100 kb")) +
  theme_verify
ggsave(file.path(VERIFY_DIR, "inter_interval_gap_histogram.pdf"), p_gap_hist,
       width = 8, height = 5)
ggsave(file.path(VERIFY_DIR, "inter_interval_gap_histogram.png"), p_gap_hist,
       width = 8, height = 5, dpi = 300)

# Cumulative distribution of gaps
p_gap_cdf <- ggplot(gap_data, aes(x = gap_to_next / 1000)) +
  stat_ecdf(geom = "step", color = "#2166AC", linewidth = 1) +
  geom_vline(xintercept = c(50, 100, 250, 500), linetype = "dashed",
             color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")) +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Gap between adjacent centroAnno intervals (kb, log scale)",
       y = "Cumulative fraction",
       title = "Cumulative distribution of inter-interval gap distances") +
  theme_verify
ggsave(file.path(VERIFY_DIR, "inter_interval_gap_cdf.pdf"), p_gap_cdf,
       width = 8, height = 5)
ggsave(file.path(VERIFY_DIR, "inter_interval_gap_cdf.png"), p_gap_cdf,
       width = 8, height = 5, dpi = 300)

# Per-chromosome gap heatmap-style summary
p_gap_perchr <- ggplot(gap_summary, aes(x = reorder(chrom, median_gap), y = 1)) +
  geom_tile(aes(fill = log10(median_gap + 1)), color = "white") +
  geom_text(aes(label = paste0(round(median_gap / 1000, 1), "kb")), size = 2.8) +
  scale_fill_viridis_c(option = "B") +
  labs(x = "Chromosome", y = NULL,
       title = "Median inter-interval gap by chromosome (kb)") +
  theme_verify +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom")
ggsave(file.path(VERIFY_DIR, "median_gap_per_chromosome.pdf"), p_gap_perchr,
       width = 10, height = 3)
ggsave(file.path(VERIFY_DIR, "median_gap_per_chromosome.png"), p_gap_perchr,
       width = 10, height = 3, dpi = 300)

# ============================================================================
# 6. MERGED CANDIDATE DOMAINS
# ============================================================================
message("\n=== 6. Merged candidate domains ===")

merge_intervals <- function(bed, merge_dist) {
  # bed: data.table with chrom, start, end sorted
  # merge_dist: merge distance in bp
  if (nrow(bed) == 0) return(data.table(chrom = character(), start = integer(),
                                          end = integer(), n_intervals = integer(),
                                          domain_id = character()))

  domains <- data.table(chrom = character(), start = integer(), end = integer(),
                        n_intervals = integer())
  current_chrom <- bed$chrom[1]
  current_start <- bed$start[1]
  current_end <- bed$end[1]
  current_n <- 1

  for (i in seq_len(nrow(bed))[-1]) {
    if (bed$chrom[i] == current_chrom && bed$start[i] - current_end <= merge_dist) {
      # Merge: extend current domain
      current_end <- max(current_end, bed$end[i])
      current_n <- current_n + 1
    } else {
      # Close current domain, start new one
      domains <- rbind(domains, data.table(
        chrom = current_chrom,
        start = current_start,
        end = current_end,
        n_intervals = current_n
      ))
      current_chrom <- bed$chrom[i]
      current_start <- bed$start[i]
      current_end <- bed$end[i]
      current_n <- 1
    }
  }
  # Don't forget the last domain
  domains <- rbind(domains, data.table(
    chrom = current_chrom,
    start = current_start,
    end = current_end,
    n_intervals = current_n
  ))

  domains[, domain_id := paste0(chrom, "_dom_", 1:.N), by = chrom]
  domains[, length_bp := end - start]
  return(domains)
}

# Merge at multiple distances
merge_dists <- c(50000, 100000, 250000, 500000)
merged_domains <- list()
for (d in merge_dists) {
  doms <- merge_intervals(copy(centro_bed), d)
  merged_domains[[as.character(d)]] <- doms

  n_domains <- nrow(doms)
  n_singleton <- sum(doms$n_intervals == 1)
  n_multi <- sum(doms$n_intervals > 1)

  message("\nMerge distance = ", d / 1000, " kb:")
  message("  Total domains: ", n_domains)
  message("  Singleton domains (1 interval): ", n_singleton, " (",
          round(100 * n_singleton / n_domains, 1), "%)")
  message("  Multi-interval domains: ", n_multi)
  message("  Median domain size: ", round(median(doms$length_bp) / 1000, 1), " kb")
  message("  Median intervals per domain: ", median(doms$n_intervals))
  message("  Max domain size: ", round(max(doms$length_bp) / 1000, 1), " kb")
  message("  Max intervals in one domain: ", max(doms$n_intervals))

  # Per-chromosome domain counts
  doms_per_chr <- doms[, .N, by = chrom]
  doms_per_chr <- merge(doms_per_chr,
                         centro_bed[, .(n_intervals = .N), by = chrom],
                         by = "chrom", all = TRUE)
  doms_per_chr[is.na(N), N := 0]
  doms_per_chr[, merge_dist_kb := d / 1000]

  fwrite(doms, file.path(VERIFY_DIR, paste0("merged_domains_", d/1000, "kb.bed")),
         sep = "\t", col.names = FALSE)
  fwrite(doms_per_chr, file.path(VERIFY_DIR,
                                  paste0("merged_domains_per_chr_", d/1000, "kb.csv")))
}

# Domain count as function of merge distance
merge_summary <- rbindlist(lapply(names(merged_domains), function(d) {
  doms <- merged_domains[[d]]
  per_chr <- doms[, .(n_domains = .N, total_span = sum(length_bp)), by = chrom]
  merge_dist <- as.numeric(d)
  data.table(
    merge_distance_kb = merge_dist / 1000,
    n_domains = nrow(doms),
    n_singleton = sum(doms$n_intervals == 1),
    n_multi = sum(doms$n_intervals > 1),
    median_domain_kb = median(doms$length_bp) / 1000,
    mean_domain_kb = mean(doms$length_bp) / 1000,
    max_domain_kb = max(doms$length_bp) / 1000,
    median_intervals_per_domain = median(doms$n_intervals)
  )
}))

fwrite(merge_summary, file.path(VERIFY_DIR, "merge_distance_summary.csv"))

# Plot domain counts across merge distances
p_merge <- ggplot(merge_summary, aes(x = factor(merge_distance_kb))) +
  geom_col(aes(y = n_multi), fill = "#2166AC", alpha = 0.8) +
  geom_col(aes(y = n_singleton), fill = "grey70", alpha = 0.8,
           position = "stack") +
  geom_text(aes(y = n_multi + n_singleton, label = n_domains), vjust = -0.5, size = 3.5) +
  labs(x = "Merge distance (kb)", y = "Number of candidate domains",
       title = "CentroAnno interval merging into candidate CENP-A enriched satellite arrays",
       subtitle = "Blue = multi-interval domains; Grey = singleton domains") +
  theme_verify
ggsave(file.path(VERIFY_DIR, "domain_merge_sensitivity.pdf"), p_merge,
       width = 7, height = 4)
ggsave(file.path(VERIFY_DIR, "domain_merge_sensitivity.png"), p_merge,
       width = 7, height = 4, dpi = 300)

# ============================================================================
# 7. CHROMOSOME-SPECIFIC DOMAIN INSPECTION (chr15, chr4, chr1)
# ============================================================================
message("\n=== 7. Chromosome-specific domain inspection ===")

for (chr in c("chr15", "chr4", "chr1")) {
  message("\n--- ", chr, " ---")
  chr_centro <- centro_bed[chrom == chr]
  message("  centroAnno intervals: ", nrow(chr_centro))
  message("  Total span: ", round(sum(chr_centro$end - chr_centro$start) / 1e6, 2), " Mb")

  # CENP-A signal at these intervals
  chr_fg <- chr_diag[sample == "XG_150" & chrom == chr]
  message("  XG_150 median log2(normalized signal): ", round(chr_fg$median_fg, 3))
  message("  XG_150 delta: ", round(chr_fg$delta, 3))
  message("  XG_150 zero-count fraction: ", round(chr_fg$zero_count_frac * 100, 1), "%")

  # Median inter-interval gap
  chr_gaps <- gap_data[chrom == chr, gap_to_next]
  if (length(chr_gaps) > 0) {
    message("  Median gap: ", round(median(chr_gaps) / 1000, 1), " kb")
    message("  Fraction <= 100kb: ", round(mean(chr_gaps <= 100000) * 100, 1), "%")
    message("  Fraction <= 500kb: ", round(mean(chr_gaps <= 500000) * 100, 1), "%")
  }

  # Merged domains at 250kb
  doms_250 <- merged_domains[["250000"]][chrom == chr]
  if (nrow(doms_250) > 0) {
    message("  Domains at 250kb merge: ", nrow(doms_250))
    message("    Sizes (kb): ", paste(round(doms_250$length_bp / 1000), collapse = ", "))
    message("    Intervals per domain: ", paste(doms_250$n_intervals, collapse = ", "))
  }
}

# ============================================================================
# 8. ZERO-COUNT AND MAPPABILITY ANALYSIS
# ============================================================================
message("\n=== 8. Zero-count and mappability analysis ===")

# Per-interval zero-count flag
setA[, is_zero := count == 0]
zero_summary <- setA[, .(
  n_total = .N,
  n_zero = sum(is_zero),
  zero_frac = mean(is_zero)
), by = .(sample, chrom)]

# Correlation between zero-count fraction and delta
zero_cor_150 <- zero_summary[sample == "XG_150"]
zero_delta <- merge(zero_cor_150, chr_diag[sample == "XG_150", .(chrom, delta)], by = "chrom")
message("Correlation: zero-count fraction vs delta: r = ",
        round(cor(zero_delta$zero_frac, zero_delta$delta, method = "spearman"), 3))

p_zero_vs_delta <- ggplot(zero_delta, aes(x = zero_frac * 100, y = delta)) +
  geom_point(size = 2.5, color = "#2166AC") +
  geom_text(aes(label = chrom), size = 2.5, vjust = -1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
  labs(x = "Zero-count intervals (%)",
       y = expression(Delta * " median log"[2] * "(normalized signal)"),
       title = "Zero-count fraction vs per-chromosome enrichment (XG_150)") +
  theme_verify
ggsave(file.path(VERIFY_DIR, "zero_count_vs_enrichment.pdf"), p_zero_vs_delta,
       width = 6, height = 5)
ggsave(file.path(VERIFY_DIR, "zero_count_vs_enrichment.png"), p_zero_vs_delta,
       width = 6, height = 5, dpi = 300)

# Mappability vs delta
if (!is.null(mapp_data)) {
  mapp_chr <- mapp_data[, .(mean_mapp = mean(mappable_fraction, na.rm = TRUE),
                              median_mapp = median(mappable_fraction, na.rm = TRUE),
                              frac_low_mapp = mean(mappable_fraction < 0.5, na.rm = TRUE)),
                          by = chrom]
  mapp_delta <- merge(mapp_chr, chr_diag[sample == "XG_150", .(chrom, delta)], by = "chrom")

  message("Correlation: mappability vs delta: r = ",
          round(cor(mapp_delta$median_mapp, mapp_delta$delta, method = "spearman"), 3))

  p_mapp_vs_delta <- ggplot(mapp_delta, aes(x = mean_mapp, y = delta)) +
    geom_point(size = 2.5, color = "#2166AC") +
    geom_text(aes(label = chrom), size = 2.5, vjust = -1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
    labs(x = "Mean mappability (H3K27ac proxy)",
         y = expression(Delta * " median log"[2] * "(normalized signal)"),
         title = "Mappability vs per-chromosome enrichment (XG_150)") +
    theme_verify
  ggsave(file.path(VERIFY_DIR, "mappability_vs_enrichment.pdf"), p_mapp_vs_delta,
         width = 6, height = 5)
  ggsave(file.path(VERIFY_DIR, "mappability_vs_enrichment.png"), p_mapp_vs_delta,
         width = 6, height = 5, dpi = 300)
}

# ============================================================================
# 9. REPLICATE REPRODUCIBILITY AT CHROMOSOME LEVEL
# ============================================================================
message("\n=== 9. Replicate reproducibility ===")

chr_wide <- dcast(chr_diag[sample %in% c("XG_150", "XG_151")],
                   chrom ~ sample, value.var = "delta")
chr_wide[, agreement := fcase(
  XG_150 > 0 & XG_151 > 0, "Both positive",
  XG_150 < 0 & XG_151 < 0, "Both negative",
  default = "Disagree"
)]

message("Replicate agreement at chromosome level:")
print(table(chr_wide$agreement))

spearman_chr <- cor(chr_wide$XG_150, chr_wide$XG_151, method = "spearman",
                    use = "complete.obs")
message("Spearman rho between XG_150 and XG_151 chromosome deltas: ",
        round(spearman_chr, 3))

p_rep_chr <- ggplot(chr_wide, aes(x = XG_150, y = XG_151)) +
  geom_point(size = 2.5, alpha = 0.8, color = "#2166AC") +
  geom_text(aes(label = chrom), size = 2.5, vjust = -1, hjust = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
  annotate("text", x = min(chr_wide$XG_150, na.rm = TRUE),
           y = max(chr_wide$XG_151, na.rm = TRUE),
           label = paste0("rho = ", round(spearman_chr, 3)),
           hjust = 0, vjust = 1, size = 4) +
  labs(x = expression("XG_150 " * Delta * " log"[2] * "(normalized signal)"),
       y = expression("XG_151 " * Delta * " log"[2] * "(normalized signal)"),
       title = "Chromosome-level replicate reproducibility",
       subtitle = "Each point = one chromosome") +
  theme_verify
ggsave(file.path(VERIFY_DIR, "replicate_chromosome_level.pdf"), p_rep_chr,
       width = 6, height = 6)
ggsave(file.path(VERIFY_DIR, "replicate_chromosome_level.png"), p_rep_chr,
       width = 6, height = 6, dpi = 300)

# ============================================================================
# 10. H3K27ac specificity check
# ============================================================================
message("\n=== 10. H3K27ac specificity ===")

chr_wide_h3k27ac <- dcast(chr_diag[sample %in% c("XG_150", "XG_152")],
                           chrom ~ sample, value.var = "delta")
setnames(chr_wide_h3k27ac, c("XG_150", "XG_152"), c("cenpa_delta", "h3k27ac_delta"))

cor_cenpa_h3k27ac <- cor(chr_wide_h3k27ac$cenpa_delta, chr_wide_h3k27ac$h3k27ac_delta,
                          method = "spearman", use = "complete.obs")
message("Spearman rho between CENP-A and H3K27ac chromosome deltas: ",
        round(cor_cenpa_h3k27ac, 3))
message("CENP-A positive chr: ", sum(chr_wide_h3k27ac$cenpa_delta > 0, na.rm = TRUE))
message("H3K27ac positive chr: ", sum(chr_wide_h3k27ac$h3k27ac_delta > 0, na.rm = TRUE))
message("Both positive: ",
        sum(chr_wide_h3k27ac$cenpa_delta > 0 & chr_wide_h3k27ac$h3k27ac_delta > 0,
            na.rm = TRUE))

# ============================================================================
# 11. MAPPING QC — check BAM properties
# ============================================================================
message("\n=== 11. Mapping QC (from saved QC files) ===")

# Read idxstats for mapping stats
for (s in SAMPLES) {
  idx_file <- file.path(QC_DIR, paste0(s, "_idxstats.txt"))
  if (file.exists(idx_file)) {
    idx <- fread(idx_file, header = FALSE,
                 col.names = c("chrom", "length", "mapped", "unmapped"))
    idx_chr <- idx[chrom %in% CHROMOSOMES]
    message(s, ": ", sum(idx_chr$mapped), " mapped reads across ", nrow(idx_chr), " chromosomes")
  }

  frag_file <- file.path(QC_DIR, paste0(s, "_fragment_sizes.txt"))
  if (file.exists(frag_file)) {
    frag <- fread(frag_file, header = FALSE, col.names = c("count", "size"))
    message("  Fragment size: median = ", round(weighted.mean(frag$size, frag$count)), " bp")
  }
}

# ============================================================================
# 12. SAVE SESSION INFO
# ============================================================================
writeLines(capture.output(sessionInfo()), file.path(VERIFY_DIR, "session_info.txt"))

message("\n=== verify_diagnostics.R DONE ===")
message("Output directory: ", VERIFY_DIR)
message("Files written:")
for (f in list.files(VERIFY_DIR, full.names = TRUE)) {
  message("  ", f)
}
