#!/usr/bin/env Rscript
# ============================================================================
# 9_methylation_plot.R -- CpG methylation level around the peak-based CENP-A
# cores (WGBS degu_6834_PFC_1, CGN-both), in the repeat-occupancy intervals.
#
# Main panel: per-chromosome mean methylation fraction by spatial interval,
# ordered short-distal -> core -> long-distal -> background; points = per
# chromosome, box = distribution, diamonds = genome-wide pooled mean.
# Supp panel: coverage (mean_cov, log10) by interval, showing the sparse but
# high-depth signature of the satellite-junction core.
#
# Usage: Rscript 9_methylation_plot.R <work_dir>
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
RES  <- file.path(WORK, "results")
MAIN <- file.path(WORK, "plots", "main")
SUPP <- file.path(WORK, "plots", "supp")
dir.create(MAIN, showWarnings = FALSE); dir.create(SUPP, showWarnings = FALSE)

m <- fread(file.path(RES, "methylation_around_cenpa_core.csv"))

REGION_ORDER <- c("short_distal", "short_intermediate", "short_proximal", "core",
                  "long_proximal", "long_intermediate", "long_distal", "background")
m[, region_key := ifelse(side %in% c("core", "background"), side,
                         paste(side, region, sep = "_"))]
m[, region_key := factor(region_key, levels = REGION_ORDER)]
REGION_LAB <- c(
  "short_distal"        = "short\ndistal 2-5 Mb",
  "short_intermediate"  = "short\nintermediate 500 kb-2 Mb",
  "short_proximal"      = "short\nproximal 0-500 kb",
  "core"                = "CENP-A\ncore",
  "long_proximal"       = "long\nproximal 0-500 kb",
  "long_intermediate"   = "long\nintermediate 500 kb-2 Mb",
  "long_distal"         = "long\ndistal 2-5 Mb",
  "background"          = "background")
m[, region_lab := REGION_LAB[as.character(region_key)]]
m[, region_lab := factor(region_lab, levels = REGION_LAB[REGION_ORDER])]

REGION_COLS <- c(short_distal = "#a6bddb", short_intermediate = "#74a9cf",
                 short_proximal = "#2b8cbe", core = "#d7191c",
                 long_proximal = "#fdae61", long_intermediate = "#f46d43",
                 long_distal = "#d73027", background = "#bdbdbd")

per_chr <- m[chromosome != "genome"]
gw      <- m[chromosome == "genome"]

# how many chromosomes have usable (>=10 cov>=5 CpGs) data in the core?
core_n <- per_chr[region_key == "core"]
cat(sprintf("Core: %d/%d chromosomes with >=10 covered CpGs; %d with 0 covered CpGs\n",
            sum(core_n$n_cpgs_cov5 >= 10), nrow(core_n),
            sum(core_n$n_cpgs_cov5 == 0)))

# ---- Panel A: methylation fraction by interval ----
pA <- ggplot(per_chr, aes(x = region_lab, y = mean_frac * 100)) +
  geom_boxplot(fill = "grey92", outlier.shape = NA, width = 0.55) +
  geom_jitter(aes(color = region_lab), width = 0.14, size = 0.9, alpha = 0.45) +
  geom_point(data = gw, aes(x = region_lab, y = mean_frac * 100),
             shape = 23, size = 2.4, fill = "black", color = "white",
             inherit.aes = FALSE) +
  scale_color_manual(values = REGION_COLS, guide = "none") +
  scale_x_discrete(labels = levels(per_chr$region_lab)) +
  labs(title = "CpG methylation around the CENP-A core (peak-based, per chromosome)",
       subtitle = "WGBS degu_6834 PFC, CpG both strands; diamonds = genome-wide pooled mean",
       x = "interval from CENP-A core edge (short/long arm)",
       y = "mean methylation fraction (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5))
ggsave(file.path(MAIN, "panel_methylation_around_cenpa_core.pdf"), pA, width = 7, height = 4.5)
ggsave(file.path(MAIN, "panel_methylation_around_cenpa_core.png"), pA, width = 7, height = 4.5, dpi = 300)

# ---- Supp panel: coverage by interval (core sparsity context) ----
pS <- ggplot(per_chr, aes(x = region_lab, y = mean_cov)) +
  geom_boxplot(fill = "grey92", outlier.shape = NA, width = 0.55) +
  geom_jitter(color = "grey40", width = 0.14, size = 0.9, alpha = 0.5) +
  geom_point(data = gw, aes(x = region_lab, y = mean_cov),
             shape = 23, size = 2.4, fill = "black", color = "white",
             inherit.aes = FALSE) +
  scale_y_log10() +
  labs(title = "WGBS coverage by interval (log10) - core is sparse but high-depth",
       subtitle = "points = per chromosome; diamonds = genome-wide pooled",
       x = "interval from CENP-A core edge", y = "mean per-CpG coverage (log10)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5))
ggsave(file.path(SUPP, "supp_methylation_coverage_by_interval.pdf"), pS, width = 7, height = 4.5)

cat("DONE 9_methylation_plot.R\n")
