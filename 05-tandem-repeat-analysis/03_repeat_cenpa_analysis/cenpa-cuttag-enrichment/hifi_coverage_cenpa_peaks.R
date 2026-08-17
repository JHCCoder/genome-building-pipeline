#!/usr/bin/env Rscript
# HiFi coverage: genome-wide average vs CENP-A peak regions vs seq31
# Evidence for collapsed/under-assembled centromere sequence.
#
# Data source: mosdepth per-base primary-reads HiFi coverage
#   hifiasm_041425_scaffolded_juiceBox_sorted_primaryReads_long_read_coverage_basepair_level.per-base.bed.gz
#   (chr4 length 151,562,733 = matches the CENP-A mapping assembly)
# Peak windows = top-1 CENP-A windows from results/top_window_signal_fraction.csv (XG_150).
# seq31 peak = the 3.9k x locus at seq31:1,218,000-1,238,000 (unplaced contig).
#
# Run: Rscript hifi_coverage_cenpa_peaks.R
# Env: r-visualizations

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# ---- Coverage values (computed via tabix on the per-base bed.gz) ----
dat <- data.frame(
  label = c("Genome-wide\naverage (30 chr)",
            "chr4\nCENP-A peak",
            "seq31\nCENP-A peak (unplaced)"),
  coverage = c(35.6, 97.6, 3920.4),
  stringsAsFactors = FALSE
)

# Order: genome avg first, then peaks ascending
dat$label <- factor(dat$label, levels = dat$label)
dat$hjust <- c(0.5, 0.5, 0.5)

p <- ggplot(dat, aes(x = label, y = coverage, fill = label)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(coverage, 1), "×")),
            vjust = -0.4, size = 4.2, fontface = "bold") +
  scale_y_log10(breaks = c(10, 30, 100, 300, 1000, 3000),
                labels = label_number(accuracy = 1)) +
  scale_fill_manual(values = c("grey55", "#4575B4", "#B2182B")) +
  labs(title = "HiFi coverage: CENP-A peaks vs genome average",
       subtitle = "Collapsed/under-assembled centromere sequence (primary-reads HiFi)",
       x = NULL,
       y = expression("mean HiFi coverage (×)"),
       caption = paste0("Genome-wide mean from mosdepth summary (30 placed chromosomes); ",
                        "peak windows = top-1 CENP-A windows; seq31 is an unplaced contig.")) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.position = "none",
    axis.text.x = element_text(size = 11, lineheight = 0.9),
    axis.title.y = element_text(size = 12),
    plot.margin = margin(10, 12, 10, 10)
  )

outbase <- "hifi_coverage_cenpa_peaks"
ggsave(paste0(outbase, ".pdf"), p, width = 6, height = 5)
ggsave(paste0(outbase, ".png"), p, width = 6, height = 5, dpi = 300)
ggsave(paste0(outbase, ".svg"), p, width = 6, height = 5, device = svglite::svglite)

cat("Wrote", paste0(outbase, ".{pdf,png,svg}"), "\n")
cat("Coverage values:\n")
print(dat[, c("label", "coverage")])
