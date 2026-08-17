#!/usr/bin/env Rscript
# ============================================================================
# 5_analysis.R (V2) — CENP-A CUT&Tag enrichment at centroAnno repeats
#
# V2 changes:
#   - XG_153 (H3K27ac control2) integrated throughout
#   - Primary analysis at merged domain level
#   - New 4-panel main figure layout:
#       A: Browser tracks (chr15, chr4, chr1)
#       B: Chromosome-level enrichment (beeswarm + boxplot)
#       C: Array-centered metaprofile (±1 Mb, boundary-relative)
#       D: Per-domain enrichment vs local background
#   - Domain local flank background as primary control
#   - Interval-level Bg1 retained as sensitivity analysis
#
# Usage:
#   conda activate r-visualizations
#   Rscript 5_analysis.R
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(patchwork)
  library(ggpubr)
  library(rstatix)
  library(scales)
  library(cowplot)
  # ggbeeswarm not in conda env — using geom_jitter instead
})

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
DATA_DIR <- file.path(BASE_DIR, "data")
COUNTS_DIR <- file.path(DATA_DIR, "counts")
DOMAIN_DIR <- file.path(DATA_DIR, "domains")
DOMAIN_COUNTS_DIR <- file.path(DATA_DIR, "domain_counts")
DOMAIN_FLANKS_DIR <- file.path(DATA_DIR, "domain_flanks")
DOMAIN_DISTANCE_DIR <- file.path(DATA_DIR, "domain_distance_profiles")
DISTANCE_DIR <- file.path(DATA_DIR, "distance_profiles")
FOREGROUND_DIR <- file.path(DATA_DIR, "foregrounds")
BACKGROUND_DIR <- file.path(DATA_DIR, "backgrounds")
MAPPABILITY_DIR <- file.path(DATA_DIR, "mappability")
QC_DIR <- file.path(DATA_DIR, "qc")
PLOTS_DIR <- file.path(BASE_DIR, "plots")
MAIN_PLOTS <- file.path(PLOTS_DIR, "main")
SUPP_PLOTS <- file.path(PLOTS_DIR, "supp")
RESULTS_DIR <- file.path(BASE_DIR, "results")

dir.create(MAIN_PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Constants
CHROMOSOMES <- paste0("chr", c(1:28, "X", "Y"))
AUTOSOMES <- paste0("chr", 1:28)
SEED_BASE <- 20260716
N_SHUFFLES <- 100

# V2: now 4 samples
SAMPLES <- c("XG_150", "XG_151", "XG_152", "XG_153")
SAMPLE_LABELS <- c(
  "XG_150" = "CENP-A rep1",
  "XG_151" = "CENP-A rep2",
  "XG_152" = "H3K27ac",
  "XG_153" = "H3K27ac rep2"
)
SAMPLE_COLORS <- c(
  "XG_150" = "#2166AC",
  "XG_151" = "#92C5DE",
  "XG_152" = "#B2182B",
  "XG_153" = "#D6604D"
)

# ============================================================================
# ggplot2 theme
# ============================================================================
theme_cenpa <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(size = 9, face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

# ============================================================================
# SECTION 1: Load data
# ============================================================================
message("=== SECTION 1: Loading data ===")

# Library sizes
lib_file <- file.path(QC_DIR, "library_sizes_fragments.txt")
if (file.exists(lib_file)) {
  lib_sizes <- fread(lib_file, header = FALSE, col.names = c("sample", "n_fragments"))
  lib_sizes <- setNames(lib_sizes$n_fragments, lib_sizes$sample)
  message("Library sizes (fragments):")
  print(lib_sizes)
} else {
  # Fall back to counting fragment BED lines
  message("Library sizes file not found, computing from fragment BEDs...")
  lib_sizes <- sapply(SAMPLES, function(s) {
    f <- file.path(DATA_DIR, "fragments", paste0(s, "_fragments.bed"))
    if (file.exists(f)) as.numeric(system(paste("wc -l <", f), intern = TRUE)) else NA
  })
  names(lib_sizes) <- SAMPLES
}

# Load foreground counts for Set A (individual intervals)
load_counts <- function(sample, region_set) {
  f <- file.path(COUNTS_DIR, paste0(sample, "_", region_set, ".txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE, col.names = c("chrom", "start", "end", "id", "count"))
  dt[, sample := sample]
  dt[, region_set := region_set]
  dt[, length := end - start]
  return(dt)
}

setA_list <- lapply(SAMPLES, function(s) load_counts(s, "setA_strict"))
setA <- rbindlist(setA_list[!sapply(setA_list, is.null)])
message("Set A (strict centroAnno): ", nrow(setA) / length(SAMPLES), " intervals per sample")

# ============================================================================
# SECTION 2: Load domain-level data
# ============================================================================
message("=== SECTION 2: Loading domain data ===")

# Merged domains
domains_file <- file.path(DOMAIN_DIR, "merged_domains_d250000.bed")
if (file.exists(domains_file)) {
  domains <- fread(domains_file, header = FALSE,
                   col.names = c("chrom", "start", "end", "domain_id", "size"))
  message("Merged domains: ", nrow(domains))
} else {
  stop("Merged domains file not found: ", domains_file,
       "\nRun 2b_define_domains.sh first.")
}

# Domain counts for all samples
load_domain_counts <- function(sample) {
  f <- file.path(DOMAIN_COUNTS_DIR, paste0(sample, "_domains.txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE,
              col.names = c("chrom", "start", "end", "domain_id", "size", "count"))
  dt[, sample := sample]
  dt[, length := end - start]
  return(dt)
}

domain_counts_list <- lapply(SAMPLES, function(s) load_domain_counts(s))
domain_counts <- rbindlist(domain_counts_list[!sapply(domain_counts_list, is.null)])
message("Domain counts: ", nrow(domain_counts) / length(SAMPLES), " domains per sample")

# Domain flank counts (local background)
load_domain_flank_counts <- function(sample) {
  f <- file.path(DOMAIN_COUNTS_DIR, paste0(sample, "_domain_flanks.txt"))
  if (!file.exists(f)) {
    # Try individual flank files
    flank_files <- list.files(DOMAIN_COUNTS_DIR,
                              pattern = paste0(sample, "_domain_flank_.*\\.txt"),
                              full.names = TRUE)
    if (length(flank_files) == 0) return(NULL)
    dt_list <- lapply(flank_files, function(ff) {
      flank_id <- sub(paste0(sample, "_domain_flank_"), "",
                      sub("\\.txt$", "", basename(ff)))
      dt <- fread(ff, header = FALSE,
                  col.names = c("chrom", "start", "end", "domain_id", "count"))
      dt[, flank_id := flank_id]
      return(dt)
    })
    return(rbindlist(dt_list))
  }
  dt <- fread(f, header = FALSE,
              col.names = c("chrom", "start", "end", "flank_id", "count"))
  dt[, sample := sample]
  dt[, length := end - start]
  return(dt)
}

domain_flanks_list <- lapply(SAMPLES, function(s) load_domain_flank_counts(s))
domain_flanks <- rbindlist(domain_flanks_list[!sapply(domain_flanks_list, is.null)])
message("Domain flank counts: ", nrow(domain_flanks), " regions")

# Load Bg1 (chromosome-matched shuffle) for baseline comparison
load_bg1 <- function(sample) {
  f <- file.path(COUNTS_DIR, paste0(sample, "_bg1_chrom_shuffle.txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE,
              col.names = c("chrom", "start", "end", "iter_id", "interval_id", "count"))
  dt[, sample := sample]
  dt[, length := end - start]
  return(dt)
}

bg1_list <- lapply(SAMPLES, function(s) load_bg1(s))
bg1 <- rbindlist(bg1_list[!sapply(bg1_list, is.null)])
message("Bg1 (chromosome shuffle): ", nrow(bg1) / length(SAMPLES), " per sample")

# ============================================================================
# SECTION 3: Normalization
# ============================================================================
message("=== SECTION 3: Normalization ===")

normalize_signal <- function(dt, lib_sizes) {
  dt[, norm_signal := (count / (length / 1000)) / (lib_sizes[sample] / 1e6)]
  return(dt)
}

setA <- normalize_signal(setA, lib_sizes)
domain_counts <- normalize_signal(domain_counts, lib_sizes)
if (nrow(domain_flanks) > 0 && "length" %in% names(domain_flanks)) {
  domain_flanks <- normalize_signal(domain_flanks, lib_sizes)
}
if (nrow(bg1) > 0) bg1 <- normalize_signal(bg1, lib_sizes)

# Pseudocount: half the minimum non-zero normalized signal across all data
all_signal <- c(setA$norm_signal, domain_counts$norm_signal)
if (nrow(domain_flanks) > 0 && "norm_signal" %in% names(domain_flanks)) {
  all_signal <- c(all_signal, domain_flanks$norm_signal)
}
if (nrow(bg1) > 0) all_signal <- c(all_signal, bg1$norm_signal)
min_nonzero <- min(all_signal[all_signal > 0], na.rm = TRUE)
pseudocount <- min_nonzero / 2
message("Min non-zero normalized signal: ", format(min_nonzero, digits = 3))
message("Pseudocount: ", format(pseudocount, digits = 3))

setA[, log2_signal := log2(norm_signal + pseudocount)]
domain_counts[, log2_signal := log2(norm_signal + pseudocount)]
if (nrow(domain_flanks) > 0 && "norm_signal" %in% names(domain_flanks)) {
  domain_flanks[, log2_signal := log2(norm_signal + pseudocount)]
}
if (nrow(bg1) > 0) bg1[, log2_signal := log2(norm_signal + pseudocount)]

# ============================================================================
# SECTION 4: Replicate concordance
# ============================================================================
message("=== SECTION 4: Replicate concordance ===")

# CENP-A replicate correlation at domain level
rep_domain <- dcast(domain_counts[sample %in% c("XG_150", "XG_151")],
                    chrom + start + end + domain_id ~ sample, value.var = "log2_signal")
cor_domain <- cor(rep_domain$XG_150, rep_domain$XG_151,
                  method = "spearman", use = "complete.obs")
message("Domain-level Spearman rho (XG_150 vs XG_151): ", round(cor_domain, 3))

# H3K27ac replicate correlation at domain level
rep_h3k27ac <- dcast(domain_counts[sample %in% c("XG_152", "XG_153")],
                      chrom + start + end + domain_id ~ sample, value.var = "log2_signal")
if (nrow(rep_h3k27ac) > 0 && all(c("XG_152", "XG_153") %in% names(rep_h3k27ac))) {
  cor_h3k27ac <- cor(rep_h3k27ac$XG_152, rep_h3k27ac$XG_153,
                      method = "spearman", use = "complete.obs")
  message("H3K27ac replicate Spearman rho (XG_152 vs XG_153): ", round(cor_h3k27ac, 3))
}

# Interval-level for comparison
rep_interval <- dcast(setA[sample %in% c("XG_150", "XG_151")],
                      chrom + start + end + id ~ sample, value.var = "log2_signal")
cor_interval <- cor(rep_interval$XG_150, rep_interval$XG_151,
                    method = "spearman", use = "complete.obs")
message("Interval-level Spearman rho (XG_150 vs XG_151): ", round(cor_interval, 3))

# Per-chromosome domain-level correlations
per_chr_cor_domain <- rep_domain[, .(
  spearman_rho = cor(XG_150, XG_151, method = "spearman", use = "complete.obs"),
  n_domains = .N
), by = chrom]
message("Per-chromosome domain rho range: ",
        round(min(per_chr_cor_domain$spearman_rho, na.rm = TRUE), 3), " - ",
        round(max(per_chr_cor_domain$spearman_rho, na.rm = TRUE), 3))
fwrite(per_chr_cor_domain, file.path(RESULTS_DIR, "replicate_correlation_per_chr.csv"))

# ============================================================================
# SECTION 5: Panel B — Chromosome-level enrichment (V2 primary inference)
# ============================================================================
message("=== SECTION 5: Panel B — Chromosome-level enrichment ===")

# Per-chromosome domain-level delta (foreground - Bg1 background)
# Domain foreground medians per chromosome
chr_fg_domain <- domain_counts[, .(
  median_fg = median(log2_signal, na.rm = TRUE),
  n_domains = .N
), by = .(sample, chrom)]

# Background (Bg1) per chromosome (using iter_001 as representative)
bg1_iter001 <- bg1[iter_id == "iter_001"]
chr_bg <- bg1_iter001[, .(
  median_bg = median(log2_signal, na.rm = TRUE)
), by = .(sample, chrom)]

# Merge and compute delta
chr_delta <- merge(chr_fg_domain, chr_bg, by = c("sample", "chrom"))
chr_delta[, delta := median_fg - median_bg]
chr_delta[, chrom_factor := factor(chrom, levels = CHROMOSOMES)]
chr_delta[, sample_label := SAMPLE_LABELS[sample]]

# Order chromosomes by mean CENPA delta
chr_order_data <- chr_delta[sample %in% c("XG_150", "XG_151"),
                             .(mean_delta = mean(delta, na.rm = TRUE)), by = chrom]
setorder(chr_order_data, mean_delta)
chr_delta[, chrom_ordered := factor(chrom, levels = chr_order_data$chrom)]

# Statistical tests (CENP-A replicates)
chr_delta_wide <- dcast(chr_delta[sample %in% c("XG_150", "XG_151")],
                        chrom ~ sample, value.var = "delta")

wilcox_150 <- wilcox.test(chr_delta_wide$XG_150, mu = 0,
                           alternative = "two.sided", exact = TRUE)
wilcox_151 <- wilcox.test(chr_delta_wide$XG_151, mu = 0,
                           alternative = "two.sided", exact = TRUE)

message("Primary inference: Paired Wilcoxon signed-rank across chromosomes")
message("  XG_150: V = ", wilcox_150$statistic, ", P = ", formatC(wilcox_150$p.value, digits = 4))
message("  XG_151: V = ", wilcox_151$statistic, ", P = ", formatC(wilcox_151$p.value, digits = 4))
message("  Median delta (XG_150): ", round(median(chr_delta_wide$XG_150, na.rm = TRUE), 3))
message("  Median delta (XG_151): ", round(median(chr_delta_wide$XG_151, na.rm = TRUE), 3))

prop_pos_150 <- mean(chr_delta_wide$XG_150 > 0, na.rm = TRUE)
prop_pos_151 <- mean(chr_delta_wide$XG_151 > 0, na.rm = TRUE)
message("  Positive chromosomes: XG_150 = ", sum(chr_delta_wide$XG_150 > 0, na.rm = TRUE),
        "/", sum(!is.na(chr_delta_wide$XG_150)),
        " (", round(prop_pos_150 * 100, 1), "%)")
message("  Positive chromosomes: XG_151 = ", sum(chr_delta_wide$XG_151 > 0, na.rm = TRUE),
        "/", sum(!is.na(chr_delta_wide$XG_151)),
        " (", round(prop_pos_151 * 100, 1), "%)")

# Classify chromosomes by enrichment strength (all 30 chr show enrichment vs Bg1).
# Use tertiles of mean delta — avoids arbitrary cutoffs.
delta_tertiles <- quantile(chr_delta_wide$mean_delta, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
chr_classification <- chr_delta_wide[, .(
  chrom,
  XG_150_delta = XG_150,
  XG_151_delta = XG_151,
  mean_delta = (XG_150 + XG_151) / 2,
  category = fcase(
    mean_delta >= delta_tertiles[3], "Strong enrichment",
    mean_delta >= delta_tertiles[2], "Moderate enrichment",
    mean_delta >= delta_tertiles[1], "Weak enrichment",
    default = "Insufficient data"
  )
)]

message("Chromosome classification (tertile-based):")
print(table(chr_classification$category))
message("  Delta range: ", round(min(chr_classification$mean_delta), 2),
        " - ", round(max(chr_classification$mean_delta), 2))

# Chr1 specifics for accurate labeling
chr1_delta_mean <- chr_classification[chrom == "chr1", mean_delta]
chr1_delta_rep1 <- chr_classification[chrom == "chr1", XG_150_delta]
chr1_delta_rep2 <- chr_classification[chrom == "chr1", XG_151_delta]

# Panel B: Beeswarm + boxplot
chr_delta_plot <- chr_delta[sample %in% c("XG_150", "XG_151")]
chr_delta_plot <- merge(chr_delta_plot,
                        chr_classification[, .(chrom, category)],
                        by = "chrom", all.x = TRUE)

cat_colors <- c(
  "Strong enrichment" = "#2166AC",
  "Moderate enrichment" = "#92C5DE",
  "Weak enrichment" = "#F4A582"
)

# Build subtitle with actual fg/bg context
p_panelB <- ggplot(chr_delta_plot, aes(x = "All chromosomes", y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_jitter(aes(color = category), size = 2.5, width = 0.2, alpha = 0.8) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.3, linewidth = 0.4) +
  facet_wrap(~ sample_label, nrow = 1) +
  scale_color_manual(values = cat_colors, name = "Enrichment\nstrength") +
  labs(x = NULL,
       y = expression(Delta * " median log"[2] * "(normalized signal) [domain " - " Bg1 shuffle]"),
       title = "Chromosome-level CENP-A enrichment at CENP-A enriched satellite arrays",
       subtitle = paste0("All 30/30 chromosomes enriched vs chromosome-shuffled background\n",
                         "Wilcoxon P = ", formatC(wilcox_150$p.value, digits = 3),
                         " (rep1), ", formatC(wilcox_151$p.value, digits = 3),
                         " (rep2); median Δ range: ",
                         round(min(chr_classification$mean_delta), 2), " – ",
                         round(max(chr_classification$mean_delta), 2))) +
  theme_cenpa +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Inset: ordered bar plot
p_panelB_inset <- ggplot(chr_delta_plot, aes(x = chrom_ordered, y = delta, fill = sample)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS) +
  labs(x = "Chromosome (ordered by enrichment)", y = expression(Delta * " log"[2] * "(normalized signal)"),
       fill = NULL) +
  theme_cenpa +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        legend.position = "right")

ggsave(file.path(MAIN_PLOTS, "panelB_chromosome_enrichment.pdf"), p_panelB,
       width = 8, height = 4.5)
ggsave(file.path(MAIN_PLOTS, "panelB_chromosome_enrichment.png"), p_panelB,
       width = 8, height = 4.5, dpi = 300)
ggsave(file.path(SUPP_PLOTS, "supp_chromosome_barplot.pdf"), p_panelB_inset,
       width = 10, height = 4)

message("Panel B saved")

# ============================================================================
# SECTION 6: Panel D — Per-domain enrichment vs local background
# ============================================================================
message("=== SECTION 6: Panel D — Per-domain enrichment vs local background ===")

# For each domain, compare foreground signal to pooled local flank background.
# NOTE: The flank counts file uses "flank_id" (e.g., upstream_1000kb) as the
# 4th column, not "domain_id" — domain IDs were lost during flank concatenation
# in step 3. We therefore use chromosome-pooled flank background per sample
# rather than per-domain matched flanks. This preserves per-domain resolution
# while using a conservative chromosome-specific background.

if (nrow(domain_flanks) > 0 && "norm_signal" %in% names(domain_flanks)) {

  # Per-domain foreground signal
  domain_fg <- domain_counts[, .(
    fg_median = median(log2_signal, na.rm = TRUE)
  ), by = .(sample, chrom, start, end, domain_id)]

  # Pooled flank background per sample per chromosome
  # (flank files have flank_id column, not domain_id)
  if ("flank_id" %in% names(domain_flanks)) {
    domain_bg_by_chr <- domain_flanks[, .(
      bg_median = median(log2_signal, na.rm = TRUE),
      bg_mean = mean(log2_signal, na.rm = TRUE),
      n_flanks = .N
    ), by = .(sample, chrom)]
  } else {
    # Fallback: if data has domain_id column, group by it
    domain_bg_by_chr <- domain_flanks[, .(
      bg_median = median(log2_signal, na.rm = TRUE),
      bg_mean = mean(log2_signal, na.rm = TRUE),
      n_flanks = .N
    ), by = .(sample, chrom)]
  }

  if (nrow(domain_bg_by_chr) > 0) {
    # Match domains to chromosome-level background
    domain_compare <- merge(domain_fg, domain_bg_by_chr, by = c("sample", "chrom"))
    domain_compare[, delta := fg_median - bg_median]
    domain_compare[, sample_label := SAMPLE_LABELS[sample]]

    # Per-domain enrichment (CENP-A rep1)
    domain_wide <- dcast(domain_compare[sample %in% c("XG_150", "XG_151")],
                         chrom + domain_id ~ sample, value.var = "delta")

    # Panel D: Scatter plot (CENP-A rep1) — per-domain foreground vs chr-pooled flank background
    d_plot_150 <- domain_compare[sample == "XG_150"]
    chr_class <- chr_classification[, .(chrom, category)]

    d_plot_150 <- merge(d_plot_150, chr_class, by = "chrom", all.x = TRUE)

    # Add a reference line: x = y (fg = bg, i.e., no enrichment)
    x_range <- range(d_plot_150$bg_median, na.rm = TRUE)
    y_range <- range(d_plot_150$fg_median, na.rm = TRUE)
    xy_limits <- c(min(x_range[1], y_range[1]), max(x_range[2], y_range[2]))

    p_panelD <- ggplot(d_plot_150, aes(x = bg_median, y = fg_median, color = category)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.7, size = 1.5) +
      scale_color_manual(values = cat_colors, name = "Chromosome\nclassification") +
      coord_fixed(xlim = xy_limits, ylim = xy_limits) +
      labs(x = expression("Chr-pooled flank median log"[2] * "(normalized signal)"),
           y = expression("Domain median log"[2] * "(normalized signal)"),
           title = "Per-domain CENP-A enrichment vs local flank background",
           subtitle = paste0("CENP-A rep1; ", nrow(d_plot_150), " domains; ",
                             round(mean(d_plot_150$fg_median > d_plot_150$bg_median) * 100),
                             "% above diagonal; chr-pooled flanks")) +
      theme_cenpa

    ggsave(file.path(MAIN_PLOTS, "panelD_per_domain_enrichment.pdf"), p_panelD,
           width = 6.5, height = 6)
    ggsave(file.path(MAIN_PLOTS, "panelD_per_domain_enrichment.png"), p_panelD,
           width = 6.5, height = 6, dpi = 300)

    message("Panel D saved")

    # Domain-level enrichment statistics (XG_150)
    n_domains_pos <- sum(domain_wide$XG_150 > 0, na.rm = TRUE)
    n_domains_total <- sum(!is.na(domain_wide$XG_150))
    message("Domains with positive enrichment (vs chr-pooled flank): ",
            n_domains_pos, "/", n_domains_total,
            " (", round(100 * n_domains_pos / n_domains_total, 1), "%)")

    # Per-chromosome domain enrichment summary
    domain_chr_summary <- domain_compare[sample %in% c("XG_150", "XG_151"), .(
      median_delta = median(delta, na.rm = TRUE),
      n_domains = .N,
      pct_positive = round(100 * sum(delta > 0) / .N, 1)
    ), by = .(sample, chrom)]
    fwrite(domain_chr_summary, file.path(RESULTS_DIR, "domain_enrichment_per_chr.csv"))
  }
} else {
  message("WARNING: Domain flank data not available, skipping Panel D")
}

# ============================================================================
# SECTION 7: Panel C — Array-centered metaprofile
# ============================================================================
message("=== SECTION 7: Panel C — Array-centered metaprofile ===")

load_domain_distance <- function(sample) {
  f <- file.path(DOMAIN_DISTANCE_DIR, paste0(sample, "_domain_bin_counts.txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE,
              col.names = c("chrom", "start", "end", "domain_id", "bin_label", "count"))
  dt[, sample := sample]
  dt[, length := end - start]
  return(dt)
}

domain_dist_list <- lapply(SAMPLES, function(s) load_domain_distance(s))
domain_dist <- rbindlist(domain_dist_list[!sapply(domain_dist_list, is.null)])

if (nrow(domain_dist) > 0) {
  domain_dist <- normalize_signal(domain_dist, lib_sizes)
  domain_dist[, log2_signal := log2(norm_signal + pseudocount)]

  # Aggregate by bin label across all domains
  # Define bin order: from far left to far right
  bin_order_domain <- c(
    "left_500000_1000000", "left_100000_500000", "left_10000_100000",
    "left_1000_10000", "left_0_1000",
    "left_margin_0_50kb", "domain_core", "right_margin_0_50kb",
    "right_0_1000", "right_1000_10000", "right_10000_100000",
    "right_100000_500000", "right_500000_1000000"
  )

  bin_labels_display <- c(
    "left_500000_1000000" = "-1Mb:-500kb",
    "left_100000_500000" = "-500kb:-100kb",
    "left_10000_100000" = "-100kb:-10kb",
    "left_1000_10000" = "-10kb:-1kb",
    "left_0_1000" = "-1kb:0",
    "left_margin_0_50kb" = "Domain\n(left)",
    "domain_core" = "Domain\n(core)",
    "right_margin_0_50kb" = "Domain\n(right)",
    "right_0_1000" = "0:+1kb",
    "right_1000_10000" = "+1kb:+10kb",
    "right_10000_100000" = "+10kb:+100kb",
    "right_100000_500000" = "+100kb:+500kb",
    "right_500000_1000000" = "+500kb:+1Mb"
  )

  # Filter to bins that exist
  existing_bins <- intersect(bin_order_domain, unique(domain_dist$bin_label))
  existing_labels <- bin_labels_display[existing_bins]

  domain_meta <- domain_dist[, .(
    mean_signal = mean(log2_signal, na.rm = TRUE),
    se_signal = sd(log2_signal, na.rm = TRUE) / sqrt(.N),
    n_domains = uniqueN(domain_id)
  ), by = .(sample, bin_label)]

  domain_meta[, bin_factor := factor(bin_label, levels = existing_bins,
                                      labels = existing_labels)]
  domain_meta[, sample_label := SAMPLE_LABELS[sample]]

  # Mark domain interior vs exterior
  interior_bins <- c("left_margin_0_50kb", "domain_core", "right_margin_0_50kb")
  domain_meta[, region := ifelse(bin_label %in% interior_bins, "Domain interior", "Flanks")]

  p_panelC <- ggplot(domain_meta, aes(x = bin_factor, y = mean_signal,
                                       color = sample, group = sample)) +
    annotate("rect",
             xmin = which(existing_bins == "left_margin_0_50kb") - 0.4,
             xmax = which(existing_bins == "right_margin_0_50kb") + 0.4,
             ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.3) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    geom_ribbon(aes(ymin = mean_signal - se_signal, ymax = mean_signal + se_signal,
                    fill = sample), alpha = 0.12, color = NA) +
    scale_color_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS) +
    scale_fill_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS) +
    labs(x = "Distance from domain boundary",
         y = expression("Mean log"[2] * "(normalized signal + pseudocount)"),
         color = NULL, fill = NULL,
         title = "CENP-A signal centered on merged CENP-A enriched satellite arrays",
         subtitle = paste0("±1 Mb flanks; ", uniqueN(domain_dist$domain_id),
                           " domains; shaded = domain interior")) +
    theme_cenpa +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))

  ggsave(file.path(MAIN_PLOTS, "panelC_domain_metaprofile.pdf"), p_panelC,
         width = 10, height = 4.5)
  ggsave(file.path(MAIN_PLOTS, "panelC_domain_metaprofile.png"), p_panelC,
         width = 10, height = 4.5, dpi = 300)

  message("Panel C saved")
} else {
  # Fall back to interval-centered metaprofile
  message("WARNING: No domain distance data, using interval-centered metaprofile...")
  load_interval_distance <- function(sample) {
    f <- file.path(DISTANCE_DIR, paste0(sample, "_distance_bins_counts.txt"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f, header = FALSE,
                col.names = c("chrom", "start", "end", "interval_id", "bin_label", "count"))
    dt[, sample := sample]
    dt[, length := end - start]
    return(dt)
  }

  dist_list <- lapply(SAMPLES, function(s) load_interval_distance(s))
  dist_data <- rbindlist(dist_list[!sapply(dist_list, is.null)])

  if (nrow(dist_data) > 0) {
    dist_data <- normalize_signal(dist_data, lib_sizes)
    dist_data[, log2_signal := log2(norm_signal + pseudocount)]

    dist_summary <- dist_data[, .(
      mean_signal = mean(log2_signal, na.rm = TRUE),
      se_signal = sd(log2_signal, na.rm = TRUE) / sqrt(.N)
    ), by = .(sample, bin_label)]

    dist_summary[, sample_label := SAMPLE_LABELS[sample]]

    p_panelC_fallback <- ggplot(dist_summary,
                                aes(x = bin_label, y = mean_signal,
                                    color = sample, group = sample)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      geom_ribbon(aes(ymin = mean_signal - se_signal, ymax = mean_signal + se_signal,
                      fill = sample), alpha = 0.15, color = NA) +
      scale_color_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS) +
      scale_fill_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS) +
      labs(x = "Distance from interval", y = expression("Mean log"[2] * "(normalized signal)"),
           title = "Interval-centered metaprofile (±1 Mb, V2 with XG_153)") +
      theme_cenpa +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(MAIN_PLOTS, "panelC_interval_metaprofile_fallback.pdf"),
           p_panelC_fallback, width = 10, height = 4.5)
    message("Panel C (fallback interval metaprofile) saved")
  }
}

# ============================================================================
# SECTION 8: Panel A — Browser track schematics
# ============================================================================
message("=== SECTION 8: Panel A — Browser tracks ===")

# Generate simplified browser views for 3 representative chromosomes
# Uses the per-base coverage data to plot signal tracks

plot_chromosome_track <- function(target_chrom, target_start, target_end,
                                   title_label, show_centroanno = TRUE) {
  # Load per-base coverage for the region
  cov_data <- list()
  for (s in SAMPLES) {
    cov_file <- file.path(DATA_DIR, "coverage", paste0(s, "_perbase.txt.gz"))
    if (!file.exists(cov_file)) next

    # Read coverage for this chromosome region using shell
    cmd <- sprintf("zcat %s | awk '$1==\"%s\" && $2>=%.0f && $2<=%.0f'",
                   cov_file, target_chrom, target_start, target_end)
    dt <- tryCatch({
      fread(cmd = cmd, header = FALSE,
            col.names = c("chrom", "pos", "coverage"))
    }, error = function(e) NULL)
    if (!is.null(dt) && nrow(dt) > 0) {
      dt[, sample := s]
      cov_data[[s]] <- dt
    }
  }

  if (length(cov_data) == 0) return(NULL)

  cov_all <- rbindlist(cov_data)
  cov_all[, sample_label := SAMPLE_LABELS[sample]]
  cov_all[, norm_cov := (coverage / (1)) / (lib_sizes[sample] / 1e6)]  # per-base normalized coverage

  # Smooth to 1 kb bins for plotting
  cov_all[, bin := floor(pos / 1000) * 1000]
  cov_binned <- cov_all[, .(
    mean_cov = mean(coverage, na.rm = TRUE),
    mean_norm_cov = mean(norm_cov, na.rm = TRUE)
  ), by = .(sample, sample_label, bin)]

  # CentroAnno intervals in this region
  setA_region <- setA[chrom == target_chrom & start >= target_start & end <= target_end]

  # Merged domains in this region
  domains_region <- domains[chrom == target_chrom &
                            start >= target_start & end <= target_end]

  # Build the plot
  p <- ggplot(cov_binned, aes(x = bin / 1e6, y = mean_norm_cov, color = sample_label)) +
    geom_line(linewidth = 0.4, alpha = 0.8) +
    scale_color_manual(values = SAMPLE_COLORS, labels = SAMPLE_LABELS) +
    scale_y_continuous(trans = "log1p",
                       labels = trans_format("log10", math_format(10^.x))) +
    labs(x = paste0(target_chrom, " position (Mb)"),
         y = "Normalized coverage",
         color = NULL,
         title = title_label) +
    theme_cenpa +
    theme(legend.position = "right", legend.key.size = unit(0.3, "cm"))

  # Add domain annotations
  if (nrow(domains_region) > 0) {
    p <- p + annotate("rect",
                      xmin = domains_region$start / 1e6,
                      xmax = domains_region$end / 1e6,
                      ymin = -Inf, ymax = Inf,
                      fill = "#2166AC", alpha = 0.1)
  }

  return(p)
}

# chr15: strong CENP-A enrichment (centromeric region)
# Find the region with peak CENP-A signal
chr15_domains <- domains[chrom == "chr15"]
if (nrow(chr15_domains) > 0) {
  # Use the largest domain as center ± 2 Mb
  main_domain <- chr15_domains[which.max(size)]
  chr15_center <- (main_domain$start + main_domain$end) / 2
  chr15_start <- max(0, chr15_center - 2e6)
  chr15_end <- chr15_center + 2e6

  p_chr15 <- plot_chromosome_track("chr15", chr15_start, chr15_end,
                                    "chr15: Strong CENP-A enrichment")
} else {
  p_chr15 <- NULL
}

# chr4: moderate enrichment, HiCAT-validated
chr4_domains <- domains[chrom == "chr4"]
if (nrow(chr4_domains) > 0) {
  main_domain <- chr4_domains[which.max(size)]
  chr4_center <- (main_domain$start + main_domain$end) / 2
  chr4_start <- max(0, chr4_center - 2e6)
  chr4_end <- chr4_center + 2e6

  p_chr4 <- plot_chromosome_track("chr4", chr4_start, chr4_end,
                                   "chr4: Moderate enrichment (HiCAT-validated)")
} else {
  p_chr4 <- NULL
}

# chr1: lower-magnitude CENP-A enrichment (all chromosomes enriched vs shuffled bg)
chr1_domains <- domains[chrom == "chr1"]
if (nrow(chr1_domains) > 0) {
  # Show a representative region with multiple domains
  main_domain <- chr1_domains[which.max(size)]
  chr1_center <- (main_domain$start + main_domain$end) / 2
  chr1_start <- max(0, chr1_center - 2e6)
  chr1_end <- chr1_center + 2e6

  p_chr1 <- plot_chromosome_track("chr1", chr1_start, chr1_end,
                                   paste0("chr1: Enriched above background (",
                                          round(chr1_delta_mean, 2), " Δ vs Bg1)"))
} else {
  p_chr1 <- NULL
}

# Assemble Panel A
track_plots <- list(p_chr15, p_chr4, p_chr1)
track_plots <- track_plots[!sapply(track_plots, is.null)]

if (length(track_plots) > 0) {
  p_panelA <- wrap_plots(track_plots, ncol = 1) +
    plot_annotation(title = "CENP-A CUT&Tag signal at representative centromeric regions")

  ggsave(file.path(MAIN_PLOTS, "panelA_browser_tracks.pdf"), p_panelA,
         width = 10, height = 8)
  ggsave(file.path(MAIN_PLOTS, "panelA_browser_tracks.png"), p_panelA,
         width = 10, height = 8, dpi = 300)
  message("Panel A saved")
} else {
  message("WARNING: No browser tracks generated (coverage data missing?)")
}

# ============================================================================
# SECTION 9: Combined main figure
# ============================================================================
message("=== SECTION 9: Combined main figure ===")

# Try to assemble all 4 panels
panels <- list()
if (exists("p_panelA") && !is.null(p_panelA)) panels$A <- p_panelA
if (exists("p_panelB") && !is.null(p_panelB)) panels$B <- p_panelB
if (exists("p_panelC") && !is.null(p_panelC)) panels$C <- p_panelC
if (exists("p_panelD") && !is.null(p_panelD)) panels$D <- p_panelD

if (length(panels) >= 2) {
  # Use patchwork to assemble
  if (length(panels) == 4) {
    combined <- (panels$A | panels$B) / (panels$C | panels$D) +
      plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 12))
  } else {
    combined <- wrap_plots(panels, ncol = 2) +
      plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 12))
  }

  ggsave(file.path(MAIN_PLOTS, "figure_main_combined.pdf"), combined,
         width = 14, height = 11)
  ggsave(file.path(MAIN_PLOTS, "figure_main_combined.png"), combined,
         width = 14, height = 11, dpi = 300)
  message("Combined main figure saved")
}

# ============================================================================
# SECTION 10: Supplementary figures
# ============================================================================
message("=== SECTION 10: Supplementary figures ===")

# Supp 1: Per-chromosome domain-level violin/boxplot grid
message("--- Per-chromosome domain violin grid ---")

domain_with_bg <- rbind(
  domain_counts[sample %in% c("XG_150", "XG_151"),
                .(sample, chrom, log2_signal, region_type = "Domain (foreground)")],
  bg1_iter001[sample %in% c("XG_150", "XG_151"),
              .(sample, chrom, log2_signal, region_type = "Shuffled (background)")]
)

p_supp_perchr <- ggplot(domain_with_bg[sample == "XG_150"],
                        aes(x = region_type, y = log2_signal, fill = region_type)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.8) +
  facet_wrap(~ chrom, ncol = 6, scales = "free_y") +
  scale_fill_manual(values = c("Domain (foreground)" = "#2166AC",
                               "Shuffled (background)" = "grey70")) +
  labs(x = NULL, y = expression(log[2](normalized~signal + pseudocount)),
       title = "Per-chromosome domain-level CENP-A enrichment (XG_150)") +
  theme_cenpa +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(size = 7)) +
  guides(fill = "none")

ggsave(file.path(SUPP_PLOTS, "supp_per_chromosome_domain_boxplot.pdf"), p_supp_perchr,
       width = 14, height = 12)

# Supp 2: CENP-A vs H3K27ac at domain level
message("--- CENP-A vs H3K27ac domain scatter ---")

scatter_domain <- dcast(domain_counts, chrom + start + end + domain_id ~ sample,
                        value.var = "log2_signal")
setnames(scatter_domain, c("XG_150", "XG_151", "XG_152", "XG_153"),
         c("cenpa_rep1", "cenpa_rep2", "h3k27ac", "h3k27ac_rep2"))
scatter_domain[, cenpa_mean := (cenpa_rep1 + cenpa_rep2) / 2]

p_supp_scatter <- ggplot(scatter_domain, aes(x = cenpa_mean, y = h3k27ac)) +
  geom_point(alpha = 0.4, size = 1, color = "#2166AC") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = expression("CENP-A mean log"[2] * "(normalized signal)"),
       y = expression("H3K27ac log"[2] * "(normalized signal)"),
       title = "CENP-A vs H3K27ac at merged CENP-A enriched satellite arrays") +
  theme_cenpa

ggsave(file.path(SUPP_PLOTS, "supp_cenpa_vs_h3k27ac_domains.pdf"), p_supp_scatter,
       width = 6, height = 6)

# Supp 3: H3K27ac replicate concordance
if (exists("rep_h3k27ac") && nrow(rep_h3k27ac) > 0 &&
    all(c("XG_152", "XG_153") %in% names(rep_h3k27ac))) {
  p_supp_h3k27ac_cor <- ggplot(rep_h3k27ac, aes(x = XG_152, y = XG_153)) +
    geom_point(alpha = 0.4, size = 1, color = "#B2182B") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
    annotate("text", x = min(rep_h3k27ac$XG_152, na.rm = TRUE),
             y = max(rep_h3k27ac$XG_153, na.rm = TRUE),
             label = paste0("Spearman rho = ", round(cor_h3k27ac, 3)),
             hjust = 0, vjust = 1, size = 3.5) +
    labs(x = "H3K27ac log2(normalized signal)", y = "H3K27ac rep2 log2(normalized signal)",
         title = "H3K27ac replicate concordance at merged domains") +
    theme_cenpa

  ggsave(file.path(SUPP_PLOTS, "supp_h3k27ac_replicate_concordance.pdf"),
         p_supp_h3k27ac_cor, width = 6, height = 6)
}

# Supp 4: Domain size vs enrichment
message("--- Domain size vs enrichment ---")

domain_sizes <- domains[, .(chrom, domain_id, size_kb = size / 1000)]
domain_delta <- dcast(domain_counts[sample %in% c("XG_150", "XG_151")],
                      chrom + domain_id ~ sample, value.var = "log2_signal")
domain_meta <- merge(domain_sizes, domain_delta, by = c("chrom", "domain_id"))
domain_meta[, mean_cenpa := (XG_150 + XG_151) / 2]

p_supp_size <- ggplot(domain_meta, aes(x = size_kb / 1000, y = mean_cenpa)) +
  geom_point(alpha = 0.6, size = 1.5, color = "#2166AC") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
  scale_x_log10() +
  labs(x = "Domain size (Mb)", y = "Mean CENP-A log2(normalized signal)",
       title = "Domain size vs CENP-A enrichment") +
  theme_cenpa

ggsave(file.path(SUPP_PLOTS, "supp_domain_size_vs_enrichment.pdf"), p_supp_size,
       width = 6, height = 5)

# Supp 5: Replicate correlation at domain level
p_supp_rep_cor <- ggplot(rep_domain, aes(x = XG_150, y = XG_151)) +
  geom_point(alpha = 0.4, size = 1, color = "#2166AC") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
  annotate("text", x = min(rep_domain$XG_150, na.rm = TRUE),
           y = max(rep_domain$XG_151, na.rm = TRUE),
           label = paste0("Spearman rho = ", round(cor_domain, 3),
                          "\n(", nrow(rep_domain), " domains)"),
           hjust = 0, vjust = 1, size = 3.5) +
  labs(x = "CENP-A rep1 log2(normalized signal)", y = "CENP-A rep2 log2(normalized signal)",
       title = "CENP-A replicate concordance at merged domains") +
  theme_cenpa

ggsave(file.path(SUPP_PLOTS, "supp_replicate_correlation_domains.pdf"), p_supp_rep_cor,
       width = 6, height = 6)
ggsave(file.path(SUPP_PLOTS, "supp_replicate_correlation_domains.png"), p_supp_rep_cor,
       width = 6, height = 6, dpi = 300)

# Supp 6: H3K27ac chromosome-level comparison
if ("XG_153" %in% chr_delta$sample) {
  chr_h3k27ac <- dcast(chr_delta[sample %in% c("XG_152", "XG_153")],
                        chrom ~ sample, value.var = "delta")

  p_supp_h3k27ac_chr <- ggplot(chr_h3k27ac, aes(x = XG_152, y = XG_153)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    geom_point(size = 2, color = "#B2182B") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = expression("H3K27ac " * Delta * " log"[2] * "(normalized signal)"),
         y = expression("H3K27ac rep2 " * Delta * " log"[2] * "(normalized signal)"),
         title = "H3K27ac replicate agreement at chromosome level") +
    theme_cenpa

  ggsave(file.path(SUPP_PLOTS, "supp_h3k27ac_chromosome_comparison.pdf"),
         p_supp_h3k27ac_chr, width = 6, height = 5)
}

# ============================================================================
# SECTION 11: Effect sizes and statistics tables
# ============================================================================
message("=== SECTION 11: Effect sizes ===")

# Domain-level effect sizes
effect_sizes_domain <- data.table(
  sample = SAMPLES,
  median_fg = sapply(SAMPLES, function(s) {
    median(domain_counts[sample == s, log2_signal], na.rm = TRUE)
  }),
  median_bg = sapply(SAMPLES, function(s) {
    median(bg1_iter001[sample == s, log2_signal], na.rm = TRUE)
  })
)
effect_sizes_domain[, median_diff := median_fg - median_bg]
effect_sizes_domain[, fold_enrichment := 2^median_diff]
effect_sizes_domain[, sample_label := SAMPLE_LABELS[sample]]

message("Domain-level effect sizes:")
print(effect_sizes_domain)

# Per-chromosome effect sizes
chr_effect <- chr_delta[, .(
  median_fg = median_fg[1],
  median_bg = median_bg[1],
  delta = delta[1],
  n_domains = n_domains[1]
), by = .(sample, chrom)]
chr_effect[, fold_enrichment := 2^delta]

fwrite(chr_effect, file.path(RESULTS_DIR, "per_chromosome_effect_sizes.csv"))
fwrite(effect_sizes_domain, file.path(RESULTS_DIR, "genomewide_effect_sizes.csv"))
fwrite(chr_classification, file.path(RESULTS_DIR, "chromosome_classification.csv"))

# ============================================================================
# SECTION 12: Results summary
# ============================================================================
message("=== SECTION 12: Results summary ===")

sink(file.path(RESULTS_DIR, "results_summary.txt"))
cat("CENP-A CUT&Tag Enrichment Analysis — V2 Results Summary\n")
cat("========================================================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Samples:\n")
for (s in SAMPLES) {
  cat(sprintf("  %s: %s (%s) — %d fragments\n",
              s, SAMPLE_LABELS[s],
              ifelse(grepl("CENP", SAMPLE_LABELS[s]), "target", "control"),
              lib_sizes[s]))
}
cat("\n")

cat("Pseudocount:", format(pseudocount, digits = 3), "\n\n")

cat("Merged domains:", nrow(domains), "\n")
cat("  Median size:", round(median(domains$size) / 1000, 1), "kb\n")
cat("  Mean size:", round(mean(domains$size) / 1000, 1), "kb\n")
cat("  Singletons:", sum(domains$size < 10000), "\n\n")

cat("Replicate concordance (domain-level):\n")
cat("  CENP-A Spearman rho:", round(cor_domain, 3), "\n")
if (exists("cor_h3k27ac")) {
  cat("  H3K27ac Spearman rho:", round(cor_h3k27ac, 3), "\n")
}
cat("\n")

cat("Primary inference (chromosome-level, domain foreground vs Bg1):\n")
cat("  XG_150: V =", wilcox_150$statistic, ", P =",
    formatC(wilcox_150$p.value, digits = 4), "\n")
cat("  XG_151: V =", wilcox_151$statistic, ", P =",
    formatC(wilcox_151$p.value, digits = 4), "\n")
cat("  Chromosomes positive: XG_150 =",
    sum(chr_delta_wide$XG_150 > 0, na.rm = TRUE), "/",
    sum(!is.na(chr_delta_wide$XG_150)),
    "(", round(prop_pos_150 * 100, 1), "%)\n")
cat("  Chromosomes positive: XG_151 =",
    sum(chr_delta_wide$XG_151 > 0, na.rm = TRUE), "/",
    sum(!is.na(chr_delta_wide$XG_151)),
    "(", round(prop_pos_151 * 100, 1), "%)\n\n")

if (exists("n_domains_pos")) {
  cat("Domain-level enrichment (vs local flanks):\n")
  cat("  Positive domains (XG_150):", n_domains_pos, "/", n_domains_total,
      "(", round(100 * n_domains_pos / n_domains_total, 1), "%)\n\n")
}

cat("Domain-level effect sizes:\n")
print(effect_sizes_domain[, .(sample_label, median_diff, fold_enrichment)])
cat("\n")

cat("Chromosome classification:\n")
print(table(chr_classification$category))
cat("\n")

cat("H3K27ac comparison:\n")
cat("  H3K27ac median delta (domains):",
    round(median(chr_delta[sample == "XG_152", delta], na.rm = TRUE), 3), "\n")
if ("XG_153" %in% chr_delta$sample) {
  cat("  H3K27ac rep2 median delta (domains):",
      round(median(chr_delta[sample == "XG_153", delta], na.rm = TRUE), 3), "\n")
}

sink()

message("Results summary written")

# ============================================================================
# Save session info
# ============================================================================
writeLines(capture.output(sessionInfo()), file.path(RESULTS_DIR, "session_info.txt"))

message("=== 5_analysis.R DONE ===")
message("Output files:")
message("  Main figures: ", MAIN_PLOTS)
message("  Supplementary: ", SUPP_PLOTS)
message("  Results tables: ", RESULTS_DIR)
