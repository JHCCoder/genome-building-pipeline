#!/usr/bin/env Rscript
# Figure 2D: PCoA of the enriched GO BP terms (novel genes vs degu background).
#
# Reads the cache from 04_prepare_figure2DE_go_clusters.R. Every significant term is a point
# sized by -log10(adj. P). The top clusters (ranked by their most significant term) are
# coloured and named in a key, in rank order; all other clusters are grey. Figure 2E colours
# its bars with the same cluster colours. Drawn at printed size for Genome Research
# (8 pt text, no title, keys inside the panel).
#
# Run with the r-clusterProfiler env:
#   Rscript 04_plot_figure2D_go_pcoa.R [width_mm] [height_mm]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

.cfg <- function(name) {
  v <- Sys.getenv(name)
  if (nzchar(v)) return(v)
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  root <- normalizePath(if (length(arg)) dirname(sub("^--file=", "", arg[1])) else getwd())
  while (!file.exists(file.path(root, "config.sh")) && root != dirname(root)) root <- dirname(root)
  cmd <- sprintf("source \"%s\" 2>/dev/null; printf \"%%s\" \"${%s}\"",
                 file.path(root, "config.sh"), name)
  v <- system2("bash", c("-c", shQuote(cmd)), stdout = TRUE)
  if (length(v) == 0 || !nzchar(v)) stop(sprintf("ERROR: %s is not set in the environment or in config.sh", name))
  v
}

args <- commandArgs(trailingOnly = TRUE)
FONT_SIZE <- 8                                    # pt at final print size; Genome Research: 8-10 pt
PANEL_W   <- (if (length(args) >= 1) as.numeric(args[1]) else 68) / 25.4   # in; ~60% of the 113 mm beside A/C
PANEL_H   <- (if (length(args) >= 2) as.numeric(args[2]) else 95) / 25.4   # in
GREY      <- "#C8C8C8"
WRAP      <- 40                                   # characters per key line

d <- read.delim(.cfg("FIG2DE_CLUSTERS_TSV"), stringsAsFactors = FALSE)
d$top <- !is.na(d$top_rank)
key <- d %>% filter(top) %>% distinct(top_rank, parentTerm, color) %>% arrange(top_rank) %>%
  mutate(key = str_wrap(parentTerm, WRAP))
d$cluster <- factor(ifelse(d$top, str_wrap(d$parentTerm, WRAP), NA), levels = key$key)
d <- d[order(d$top, d$score), ]                   # coloured and larger points drawn last

out_dir  <- .cfg("FIG2DE_OUT_DIR")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_stem <- file.path(out_dir, sprintf("figure2D_go_pcoa_%.0fx%.0fmm", PANEL_W * 25.4, PANEL_H * 25.4))

p <- ggplot(d, aes(x = PCoA1, y = PCoA2)) +
  geom_point(data = d[!d$top, ], aes(size = score), colour = GREY, stroke = 0) +
  geom_point(data = d[d$top, ], aes(size = score, colour = cluster), alpha = 0.9, stroke = 0) +
  scale_colour_manual(name = NULL, values = setNames(key$color, key$key), drop = FALSE) +
  scale_size_continuous(name = "−log10 adj. P", range = c(0.6, 4),   # plain text: a plotmath subscript prints < 8 pt
                        breaks = c(2, 5, 10, 20)) +
  labs(x = "PCoA 1", y = "PCoA 2") +
  theme_classic(base_size = FONT_SIZE, base_family = "Nimbus Sans") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = FONT_SIZE),
        axis.line = element_line(linewidth = 0.5 / .pt * 1.5),
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.justification = "left",
        legend.title = element_text(size = FONT_SIZE), legend.text = element_text(size = FONT_SIZE, lineheight = 0.95),
        legend.key.height = unit(3.2, "mm"), legend.key.width = unit(3, "mm"),
        legend.spacing.y = unit(0.5, "mm"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -8),
        plot.margin = margin(1, 1, 1, 1, "mm")) +
  guides(colour = guide_legend(order = 1, ncol = 1, override.aes = list(size = 2.6)),
         size = guide_legend(order = 2, override.aes = list(colour = "grey40"), nrow = 1))

# SVG: measure text with Nimbus Sans but write "Helvetica"
nimbus <- "/usr/share/fonts/urw-base35/NimbusSans-"
nimbus_as_helvetica <- list(
  plain = list(alias = "Helvetica", file = paste0(nimbus, "Regular.otf")),
  bold = list(alias = "Helvetica", file = paste0(nimbus, "Bold.otf")),
  italic = list(alias = "Helvetica", file = paste0(nimbus, "Italic.otf")),
  bolditalic = list(alias = "Helvetica", file = paste0(nimbus, "BoldItalic.otf")))
svg_user_fonts <- list("sans" = nimbus_as_helvetica, "Nimbus Sans" = nimbus_as_helvetica,
                       "Helvetica" = nimbus_as_helvetica)

stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
for (ext in c("pdf", "svg", "png")) {
  f <- paste0(out_stem, ".", ext)
  if (file.exists(f)) file.rename(f, paste0(f, ".pre-", stamp))
}
ggsave(paste0(out_stem, ".pdf"), p, width = PANEL_W, height = PANEL_H, units = "in", device = cairo_pdf)
ggsave(paste0(out_stem, ".svg"), p, width = PANEL_W, height = PANEL_H, units = "in", user_fonts = svg_user_fonts)
ggsave(paste0(out_stem, ".png"), p, width = PANEL_W, height = PANEL_H, units = "in", dpi = 600)
cat(out_stem, ".{pdf,svg,png}\n", sep = "")
