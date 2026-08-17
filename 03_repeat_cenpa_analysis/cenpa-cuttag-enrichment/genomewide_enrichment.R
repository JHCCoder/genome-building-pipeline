#!/usr/bin/env Rscript
# genomewide_enrichment.R — Reproduce V1 genome-wide enrichment analysis
# Compares centroAnno intervals (setA) vs shuffled background (Bg1) for all 4 samples
# Split version: higher-enrichment chromosomes vs lower-enrichment chromosomes (all enriched)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  # cliff.delta computed manually (effsize package not in conda env)
})

# ============================================================================
# Paths
# ============================================================================
BASE_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
DATA_DIR <- file.path(BASE_DIR, "data")
COUNTS_DIR <- file.path(DATA_DIR, "counts")
RESULTS_DIR <- file.path(BASE_DIR, "results")
PLOTS_DIR <- file.path(BASE_DIR, "plots/main")

# ============================================================================
# Sample metadata
# ============================================================================
SAMPLES <- c("XG_150", "XG_151", "XG_152", "XG_153")
SAMPLE_LABELS <- c(
  "XG_150" = "CENP-A rep1",
  "XG_151" = "CENP-A rep2",
  "XG_152" = "H3K27ac ctrl1",
  "XG_153" = "H3K27ac ctrl2"
)
SAMPLE_COLORS <- c(
  "XG_150" = "#2166AC",
  "XG_151" = "#92C5DE",
  "XG_152" = "#F4A582",
  "XG_153" = "#B2182B"
)

# Chromosome classification: compute tertiles from per-chromosome effect sizes
# All 30 chromosomes show enrichment vs shuffled background; classify by magnitude
chr_effects <- fread(file.path(RESULTS_DIR, "per_chromosome_effect_sizes.csv"))
chr_delta <- chr_effects[, .(
  mean_delta = mean(delta, na.rm = TRUE)
), by = chrom]
delta_tertiles <- quantile(chr_delta$mean_delta, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
chr_delta[, category := fcase(
  mean_delta >= delta_tertiles[3], "Strong enrichment",
  mean_delta >= delta_tertiles[2], "Moderate enrichment",
  mean_delta >= delta_tertiles[1], "Lower enrichment",
  default = "Lower enrichment"
)]
chr_class <- chr_delta
# All chromosomes are enriched; split top 2/3 vs bottom 1/3 for comparison
chr_class[, enriched := fifelse(
  category %in% c("Strong enrichment", "Moderate enrichment"),
  "Higher enrichment",
  "Lower enrichment"
)]
message("Chromosome classification (tertile-based, all chr enriched vs background):")
print(table(chr_class$enriched))
print(table(chr_class$category))
print(chr_class[, .(chrom, mean_delta, category, enriched)])

# ============================================================================
# Theme
# ============================================================================
theme_cenpa <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

# ============================================================================
# Load data
# ============================================================================
message("=== Loading data ===")

# Load setA counts (centroAnno strict intervals)
load_setA <- function(sample) {
  f <- file.path(COUNTS_DIR, paste0(sample, "_setA_strict.txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE,
              col.names = c("chrom", "start", "end", "id", "count"))
  dt[, sample := sample]
  dt[, length := end - start]
  return(dt)
}

# Load Bg1 counts (chromosome-shuffled background, 100 iterations)
load_bg1 <- function(sample) {
  f <- file.path(COUNTS_DIR, paste0(sample, "_bg1_chrom_shuffle.txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, header = FALSE,
              col.names = c("chrom", "start", "end", "iter_id", "interval_id", "count"))
  dt[, sample := sample]
  dt[, length := end - start]
  return(dt)
}

# Library sizes (total fragments)
lib_sizes <- c(
  "XG_150" = 27156788,
  "XG_151" = 25305681,
  "XG_152" = 30805412,
  "XG_153" = 111250740
)

setA_list <- lapply(SAMPLES, load_setA)
setA <- rbindlist(setA_list[!sapply(setA_list, is.null)])
message("Set A (centroAnno): ", nrow(setA) / length(SAMPLES), " intervals per sample")

bg1_list <- lapply(SAMPLES, load_bg1)
bg1 <- rbindlist(bg1_list[!sapply(bg1_list, is.null)])
message("Bg1 (shuffled): ", nrow(bg1) / length(SAMPLES), " intervals per sample (100 iterations)")

# ============================================================================
# Normalize CUT&Tag signal (fragments per kb per million mapped)
# ============================================================================
message("=== Normalizing CUT&Tag signal ===")

normalize_signal <- function(dt, lib_sizes) {
  dt[, norm_signal := (count / (length / 1000)) / (lib_sizes[sample] / 1e6)]
  return(dt)
}

setA <- normalize_signal(setA, lib_sizes)
bg1 <- normalize_signal(bg1, lib_sizes)

# Pseudocount: half the minimum non-zero normalized signal
all_signal <- c(setA$norm_signal, bg1$norm_signal)
min_nonzero <- min(all_signal[all_signal > 0], na.rm = TRUE)
pseudocount <- min_nonzero / 2
message("Min non-zero normalized signal: ", format(min_nonzero, digits = 4))
message("Pseudocount: ", format(pseudocount, digits = 4))

setA[, log2_signal := log2(norm_signal + pseudocount)]
bg1[, log2_signal := log2(norm_signal + pseudocount)]

# ============================================================================
# Build comparison dataset (use bg1 iter_001 as representative)
# ============================================================================
message("=== Building comparison dataset ===")

bg1_iter001 <- bg1[iter_id == "iter_001"]

compare_data <- rbind(
  setA[, .(sample, chrom, start, end, log2_signal, region_type = "centroAnno")],
  bg1_iter001[, .(sample, chrom, start, end, log2_signal, region_type = "shuffled\nbackground")]
)
compare_data[, sample_label := factor(SAMPLE_LABELS[sample],
                                       levels = unname(SAMPLE_LABELS[SAMPLES]))]

# Merge chromosome classification
compare_data <- merge(compare_data, chr_class[, .(chrom, enriched)],
                      by = "chrom", all.x = TRUE)

message("Comparison data: ", nrow(compare_data), " rows")
message("Samples: ", paste(unique(compare_data$sample_label), collapse = ", "))

# ============================================================================
# Statistics function
# ============================================================================
compute_stats <- function(data, sample_name, group_col = NULL, group_val = NULL) {
  if (!is.null(group_col)) {
    d <- data[sample == sample_name & get(group_col) == group_val]
    label_suffix <- paste0(" (", group_val, ")")
  } else {
    d <- data[sample == sample_name]
    label_suffix <- ""
  }
  fg <- d[region_type == "centroAnno", log2_signal]
  bg <- d[region_type == "shuffled\nbackground", log2_signal]

  if (length(fg) < 5 || length(bg) < 5) {
    return(data.table(
      sample = sample_name, group = paste0(sample_name, label_suffix),
      median_fg = NA_real_, median_bg = NA_real_,
      median_delta = NA_real_, p_wilcox = NA_real_,
      cliffs_delta = NA_real_, fold_enrichment = NA_real_,
      n_fg = length(fg), n_bg = length(bg)
    ))
  }

  wt <- wilcox.test(fg, bg, alternative = "two.sided")
  # Cliff's delta: proportion of fg > bg minus proportion of fg < bg
  # Range: -1 to 1; 0 = complete overlap; ±1 = complete separation
  cliff_d <- function(x, y) {
    # For efficiency with large vectors, sample if needed
    nx <- length(x); ny <- length(y)
    if (nx * ny > 1e7) {
      # Subsample each vector so total comparisons ≤ 1e6
      target_n <- min(nx, ny, 1000)
      xi <- sample(nx, target_n)
      yi <- sample(ny, target_n)
      x <- x[xi]; y <- y[yi]
    }
    gt <- outer(x, y, function(a, b) ifelse(a > b, 1, 0))
    lt <- outer(x, y, function(a, b) ifelse(a < b, 1, 0))
    (sum(gt) - sum(lt)) / length(gt)
  }
  cd_val <- cliff_d(fg, bg)

  data.table(
    sample = sample_name,
    group = paste0(sample_name, label_suffix),
    median_fg = median(fg, na.rm = TRUE),
    median_bg = median(bg, na.rm = TRUE),
    median_delta = median(fg, na.rm = TRUE) - median(bg, na.rm = TRUE),
    p_wilcox = wt$p.value,
    cliffs_delta = cd_val,
    fold_enrichment = 2^(median(fg, na.rm = TRUE) - median(bg, na.rm = TRUE)),
    n_fg = length(fg),
    n_bg = length(bg)
  )
}

# ============================================================================
# FIGURE 1: Genome-wide enrichment — all chromosomes pooled
# ============================================================================
message("=== Figure 1: Genome-wide enrichment (all chromosomes) ===")

stats_all <- rbindlist(lapply(SAMPLES, function(s) compute_stats(compare_data, s)))
stats_all[, p_label := fifelse(
  p_wilcox < 0.001, "P < 0.001",
  paste0("P = ", formatC(p_wilcox, digits = 3))
)]

message("Genome-wide statistics:")
print(stats_all[, .(sample, median_fg, median_bg, median_delta,
                     p_label, cliffs_delta, fold_enrichment)])

p1 <- ggplot(compare_data, aes(x = region_type, y = log2_signal, fill = sample)) +
  geom_violin(alpha = 0.6, linewidth = 0.3, draw_quantiles = 0.5) +
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.3, linewidth = 0.3) +
  facet_wrap(~ sample_label, nrow = 1) +
  scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
  labs(x = NULL, y = expression(log[2](normalized~signal + pseudocount)),
       title = "CENP-A enriches at centroAnno-predicted repeats genome-wide",
       subtitle = paste0("All 3,420 centroAnno intervals vs. chromosome-matched shuffled background; ",
                         nrow(compare_data) / 2, " intervals per violin")) +
  theme_cenpa

# Add P-value annotations
stats_anno <- stats_all
stats_anno[, x_pos := 1.5]
stats_anno[, y_pos := max(compare_data$log2_signal, na.rm = TRUE) * 0.95]
stats_anno[, sample_label := factor(SAMPLE_LABELS[sample],
                                     levels = unname(SAMPLE_LABELS[SAMPLES]))]

p1 <- p1 +
  geom_text(data = stats_anno, aes(x = x_pos, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 2.8, hjust = 0.5, fontface = "italic")

ggsave(file.path(PLOTS_DIR, "genomewide_enrichment_all_samples.pdf"), p1,
       width = 12, height = 5)
ggsave(file.path(PLOTS_DIR, "genomewide_enrichment_all_samples.png"), p1,
       width = 12, height = 5, dpi = 300)
message("Figure 1 saved")

fwrite(stats_all, file.path(RESULTS_DIR, "genomewide_enrichment_stats_all.csv"))
message("Stats table saved")

# ============================================================================
# FIGURE 2: Split by enrichment magnitude (all chromosomes enriched, split by tertile)
# ============================================================================
message("=== Figure 2: Split by enrichment magnitude ===")

# Compute stats for enriched vs not-enriched subsets
stats_split <- rbindlist(lapply(SAMPLES, function(s) {
  rbind(
    compute_stats(compare_data, s, "enriched", "Higher enrichment"),
    compute_stats(compare_data, s, "enriched", "Lower enrichment")
  )
}))
stats_split[, p_label := fifelse(
  is.na(p_wilcox), "N/A",
  fifelse(p_wilcox < 0.001, "P < 0.001",
          paste0("P = ", formatC(p_wilcox, digits = 3)))
)]
stats_split[, sample_label := factor(SAMPLE_LABELS[sample],
                                      levels = unname(SAMPLE_LABELS[SAMPLES]))]

message("Split statistics:")
print(stats_split[, .(sample, group, median_fg, median_bg, median_delta,
                       p_label, cliffs_delta, fold_enrichment, n_fg)])

# Prep plot data
compare_data[, enriched_label := factor(enriched,
                                         levels = c("Higher enrichment", "Lower enrichment"))]

p2 <- ggplot(compare_data, aes(x = region_type, y = log2_signal, fill = sample)) +
  geom_violin(alpha = 0.6, linewidth = 0.3, draw_quantiles = 0.5) +
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.3, linewidth = 0.3) +
  facet_grid(enriched_label ~ sample_label) +
  scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
  labs(x = NULL, y = expression(log[2](normalized~signal + pseudocount)),
       title = "CENP-A enrichment at centroAnno repeats, split by chromosome enrichment magnitude",
       subtitle = paste0("Higher enrichment: ", sum(chr_class$enriched == "Higher enrichment"),
                         " chromosomes; Lower enrichment: ",
                         sum(chr_class$enriched == "Lower enrichment"),
                         " chromosomes (all chromosomes enriched vs background)")) +
  theme_cenpa

# Add P-values
stats_split[, x_pos := 1.5]
stats_split[, y_pos := max(compare_data$log2_signal, na.rm = TRUE) * 0.90]
stats_split[, enriched_label := factor(group,
                                        levels = c("XG_150 (Higher enrichment)",
                                                   "XG_151 (Higher enrichment)",
                                                   "XG_152 (Higher enrichment)",
                                                   "XG_153 (Higher enrichment)",
                                                   "XG_150 (Lower enrichment)",
                                                   "XG_151 (Lower enrichment)",
                                                   "XG_152 (Lower enrichment)",
                                                   "XG_153 (Lower enrichment)"))]

# Simpler annotation: add p-values as a subtitle row
# Actually, faceted annotation is complex. Let's add hlines and annotations inline.
# For clarity, just add the stats as a separate table in the plot
p2 <- p2 +
  geom_text(data = stats_split,
            aes(x = x_pos, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 2.5, hjust = 0.5, fontface = "italic")

ggsave(file.path(PLOTS_DIR, "genomewide_enrichment_by_enrichment_status.pdf"), p2,
       width = 12, height = 7)
ggsave(file.path(PLOTS_DIR, "genomewide_enrichment_by_enrichment_status.png"), p2,
       width = 12, height = 7, dpi = 300)
message("Figure 2 saved")

fwrite(stats_split, file.path(RESULTS_DIR, "genomewide_enrichment_stats_split.csv"))
message("Split stats saved")

# ============================================================================
# Summary table
# ============================================================================
message("\n=== SUMMARY ===")
message("Key question: Are centroAnno-annotated sequences enriched for CENP-A?")
message("")

for (s in c("XG_150", "XG_151")) {
  st <- stats_all[sample == s]
  message(SAMPLE_LABELS[s], ":")
  message("  Genome-wide: delta = ", round(st$median_delta, 3),
          ", fold enrichment = ", round(st$fold_enrichment, 2),
          ", Cliff's delta = ", round(st$cliffs_delta, 3),
          ", P = ", formatC(st$p_wilcox, digits = 3))

  st_enr <- stats_split[sample == s & group == paste0(s, " (Higher enrichment)")]
  message("  Higher-enrichment chromosomes: delta = ", round(st_enr$median_delta, 3),
          ", fold enrichment = ", round(st_enr$fold_enrichment, 2),
          ", n_intervals = ", st_enr$n_fg)

  st_low <- stats_split[sample == s & group == paste0(s, " (Lower enrichment)")]
  message("  Lower-enrichment chromosomes: delta = ", round(st_low$median_delta, 3),
          ", fold enrichment = ", round(st_low$fold_enrichment, 2),
          ", n_intervals = ", st_low$n_fg)
}
message("")
message("Done.")
