#!/usr/bin/env Rscript
# ============================================================================
# 9b_methylation_349bp_plot.R -- CpG methylation of the 348-349 bp satellite,
# by spatial interval around the peak-based CENP-A core.
#
# Reads the cached table from 7b_methylation_349bp_profile.py
#   (results/methylation_349bp_satellite_by_region.csv)
# and the genome-wide 349-bp occupancy from the repeat-occupancy table
#   (results/repeat_occupancy_peakbased_with_genome.csv, bin 6 = 348-349 bp)
# -- i.e. the "genomewide_occupancy_by_bin" occupancy for the 349-bp bin.
#
# Single panel: per-chromosome mean methylation fraction of satellite CpGs
# (box = distribution, jittered points = per chromosome, diamond = genome-wide
# pooled mean) across the 8 spatial intervals, with the 349-bp satellite
# occupancy (%) printed under each region label.
#
# Outputs: plots/main/panel_methylation_349bp_satellite.{pdf,png,svg}
# (SVG via svglite so it opens correctly in Adobe Illustrator).
#
# Usage: Rscript 9b_methylation_349bp_plot.R <work_dir>   (r-visualizations env)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
RES  <- file.path(WORK, "results")
MAIN <- file.path(WORK, "plots", "main")
dir.create(MAIN, showWarnings = FALSE)

# ---- single font-size knob for manual tuning ----
FONT_SIZE <- 12

m <- fread(file.path(RES, "methylation_349bp_satellite_by_region.csv"))

# ---- genome-wide 349-bp occupancy per region (bin 6) for the x-axis labels ----
occ <- fread(file.path(RES, "repeat_occupancy_peakbased_with_genome.csv"))
occ6 <- occ[repeat_bin == 6 & chromosome == "genome",
            .(side, region, region_order, occ_pct = repeat_percent)]

REGION_ORDER <- c("short_distal", "short_intermediate", "short_proximal", "core",
                  "long_proximal", "long_intermediate", "long_distal", "background")
m[, region_key := ifelse(side %in% c("core", "background"), side,
                         paste(side, region, sep = "_"))]
m[, region_key := factor(region_key, levels = REGION_ORDER)]
occ6[, region_key := ifelse(side %in% c("core", "background"), side,
                            paste(side, region, sep = "_"))]

# region base label (2 lines, as in 9_methylation_plot.R) + occupancy as 3rd line
REGION_LAB <- c(
  "short_distal"        = "distal\n2-5 Mb",
  "short_intermediate"  = "intermediate\n500 kb-2 Mb",
  "short_proximal"      = "proximal\n0-500 kb",
  "core"                = "CENP-A\ncore",
  "long_proximal"       = "proximal\n0-500 kb",
  "long_intermediate"   = "intermediate\n500 kb-2 Mb",
  "long_distal"         = "distal\n2-5 Mb",
  "background"          = "background")
lab <- vapply(REGION_ORDER, function(k) {
  pct <- occ6[region_key == k, occ_pct]
  pct_s <- if (length(pct)) sprintf("%.1f%%", pct[1]) else ""
  sprintf("%s\n%s", REGION_LAB[[k]], pct_s)
}, character(1))

# prefix "short"/"long" as a top strip so the occupancy % sits at the axis
# (short arm above, long arm above, occupancy = bottom line of each tick label)
pref <- c("short_distal" = "short", "short_intermediate" = "short",
          "short_proximal" = "short", "core" = "",
          "long_proximal" = "long", "long_intermediate" = "long",
          "long_distal" = "long", "background" = "")
lab_full <- vapply(REGION_ORDER, function(k) {
  p <- pref[[k]]
  if (nzchar(p)) sprintf("%s\n%s", p, lab[[k]]) else lab[[k]]
}, character(1))

m[, region_lab := lab_full[as.character(region_key)]]
m[, region_lab := factor(region_lab, levels = lab_full[REGION_ORDER])]

SATELLITE_COL <- "#2166ac"   # 348-349 bp bin colour (project palette)

per_chr <- m[chromosome != "genome" & !is.na(mean_frac)]
gw      <- m[chromosome == "genome"]

# how many chromosomes have usable (>= MIN_CPG) satellite CpGs in the core?
core_n <- per_chr[region_key == "core"]
cat(sprintf("Core: %d/%d chromosomes with >=10 covered satellite CpGs; n=%d total\n",
            sum(core_n$n_cpgs_cov5 >= 10), nrow(core_n), sum(core_n$n_cpgs_cov5)))

# ---- main panel: methylation of satellite CpGs by interval + occupancy ----
p <- ggplot(per_chr, aes(x = region_lab, y = mean_frac * 100)) +
  geom_boxplot(fill = "grey92", outlier.shape = NA, width = 0.55) +
  geom_jitter(color = SATELLITE_COL, width = 0.14, size = 1.1, alpha = 0.55) +
  geom_point(data = gw, aes(x = region_lab, y = mean_frac * 100),
             shape = 23, size = 2.6, fill = "black", color = "white",
             inherit.aes = FALSE) +
  scale_x_discrete(labels = levels(per_chr$region_lab)) +
  labs(title = "Methylation of 349-bp satellite CpGs around the CENP-A core",
       subtitle = "WGBS degu_6834 PFC, CpG both strands, cov>=5; diamonds = genome-wide pooled mean;\nbottom % = 349-bp satellite occupancy (genome-wide pooled)",
       x = NULL,
       y = "mean methylation fraction (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = FONT_SIZE - 4),
        axis.text.y = element_text(size = FONT_SIZE - 2),
        axis.title.y = element_text(size = FONT_SIZE),
        plot.title = element_text(size = FONT_SIZE + 2, face = "bold"),
        plot.subtitle = element_text(size = FONT_SIZE - 3, color = "grey30"),
        legend.position = "none")

out_base <- file.path(MAIN, "panel_methylation_349bp_satellite")
ggsave(paste0(out_base, ".pdf"), p, width = 7, height = 4.8)
ggsave(paste0(out_base, ".png"), p, width = 7, height = 4.8, dpi = 300, bg = "white")
ggsave(paste0(out_base, ".svg"), p, width = 7, height = 4.8, bg = "white",
       device = svglite::svglite)

cat("DONE 9b_methylation_349bp_plot.R ->", out_base, "{pdf,png,svg}\n")
