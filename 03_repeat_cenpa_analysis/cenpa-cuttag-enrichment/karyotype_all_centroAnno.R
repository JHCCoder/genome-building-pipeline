#!/usr/bin/env Rscript
# ============================================================================
# karyotype_all_centroAnno.R — Karyotype with ALL centroAnno raw intervals
#
# Two tracks:
#   1. (above ideogram) CENP-A CUT&Tag signal (25-kb bins) + triangle markers
#   2. (below ideogram) All centroAnno satellite arrays (raw Set A, n=3,420)
#
# No CENP-A filtering — this shows what centroAnno calls genome-wide.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(karyoploteR)
  library(rtracklayer)
  library(GenomicRanges)
})

CHROMOSOMES <- paste0("chr", c(1:28, "X", "Y"))
BIN_WIDTH <- 25000L

# ----------------------------------------------------------------------------
# Input files
# ----------------------------------------------------------------------------
fai_file <- paste0(
  "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/",
  "degu-genome-browser-pythonVersion/",
  "assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta.fai"
)

bw_file <- paste0(
  "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/",
  "figure/cenpa-cuttag-centromere/bw_files/XG_150.all.bw"
)

BASE_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
FOREGROUND_DIR <- file.path(BASE_DIR, "data", "foregrounds")

# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------
# Raw centroAnno intervals (all 3,420, no CENP-A filtering)
centroAnno_raw <- fread(file.path(FOREGROUND_DIR, "setA_centroAnno_strict.bed"),
                        header = FALSE,
                        col.names = c("chrom", "start", "end", "id"))
centroAnno_raw <- centroAnno_raw[chrom %in% CHROMOSOMES]
cat(sprintf("centroAnno raw satellite arrays: %d\n", nrow(centroAnno_raw)))

# ----------------------------------------------------------------------------
# Chromosome sizes
# ----------------------------------------------------------------------------
chr_sizes <- fread(fai_file, header = FALSE, select = c(1, 2),
                   col.names = c("chrom", "size"))
chr_sizes <- chr_sizes[chrom %in% CHROMOSOMES]
chr_sizes[, chromosome_order := match(chrom, CHROMOSOMES)]
setorder(chr_sizes, chromosome_order)
chr_lengths <- setNames(chr_sizes$size, chr_sizes$chrom)

# ----------------------------------------------------------------------------
# Load BigWig and bin into 25-kb windows
# ----------------------------------------------------------------------------
bw <- BigWigFile(bw_file)

all_bins <- list()
for (chrom_i in CHROMOSOMES) {
  chr_len <- chr_lengths[[chrom_i]]
  n_bins <- max(1L, as.integer(ceiling(chr_len / BIN_WIDTH)))
  summarized <- summary(bw, size = n_bins, type = "mean",
                        which = GRanges(seqnames = chrom_i,
                                        ranges = IRanges(start = 1, end = chr_len)))[[1]]
  scores <- as.numeric(summarized$score)
  scores[is.na(scores) | !is.finite(scores)] <- 0
  scores <- pmax(scores, 0)
  all_bins[[chrom_i]] <- data.table(
    chrom = chrom_i,
    start = start(summarized),
    end   = end(summarized),
    score = scores
  )
}
bins_dt <- rbindlist(all_bins)

# Per-chromosome max for separate scaling
chr_signal_max <- bins_dt[, .(chr_max = max(score, na.rm = TRUE)), by = chrom]
chr_signal_max[!is.finite(chr_max) | chr_max <= 0, chr_max := 1]
chr_max_lookup <- setNames(chr_signal_max$chr_max, chr_signal_max$chrom)

# Widen bins for display
BIN_DISPLAY_MIN <- 200000L
bins_dt[, bin_midpoint := floor((start + end) / 2)]
bins_dt[, display_width := pmax(end - start + 1L, BIN_DISPLAY_MIN)]
bins_dt[, chr_len := chr_lengths[chrom]]
bins_dt[, plot_start := pmax(1L, floor(bin_midpoint - display_width / 2))]
bins_dt[, plot_end := pmin(chr_len, ceiling(bin_midpoint + display_width / 2))]

bins_positive <- bins_dt[score > 0]
bins_gr <- toGRanges(data.frame(
  chr = bins_positive$chrom,
  start = bins_positive$plot_start,
  end = bins_positive$plot_end,
  y = bins_positive$score
))

# ----------------------------------------------------------------------------
# Identify highest-mean 100-kb window per chromosome (for triangle marker)
# ----------------------------------------------------------------------------
SUSTAINED_N_BINS <- 4L  # 4 × 25 kb = 100 kb
top_regions <- lapply(CHROMOSOMES, function(chrom_i) {
  chr_dt <- bins_dt[chrom == chrom_i][order(start)]
  if (nrow(chr_dt) < SUSTAINED_N_BINS) return(NULL)
  chr_dt[, roll_mean := frollmean(score, n = SUSTAINED_N_BINS, align = "left", fill = NA)]
  valid <- which(is.finite(chr_dt$roll_mean))
  if (length(valid) == 0) return(NULL)
  best <- valid[which.max(chr_dt$roll_mean[valid])]
  end_idx <- best + SUSTAINED_N_BINS - 1L
  sel <- chr_dt[best:end_idx]
  peak_idx <- which.max(sel$score)[1]
  data.table(
    chrom = chrom_i,
    region_start = min(sel$start), region_end = max(sel$end),
    peak_midpoint = floor((sel$start[peak_idx] + sel$end[peak_idx]) / 2),
    peak_score = sel$score[peak_idx],
    chr_max = chr_max_lookup[chrom_i]
  )
})
top_regions <- rbindlist(Filter(Negate(is.null), top_regions))
top_regions[, peak_fraction := pmin(1, pmax(0, peak_score / chr_max))]

# ----------------------------------------------------------------------------
# Prepare GRanges for centroAnno track
# ----------------------------------------------------------------------------
genome_df <- data.frame(chr = CHROMOSOMES, start = 1, end = unname(chr_lengths[CHROMOSOMES]))
custom_genome <- toGRanges(genome_df)

CENTROANNO_DISPLAY_MIN <- 200000L
centroAnno_df <- as.data.frame(centroAnno_raw)
centroAnno_df$width <- centroAnno_df$end - centroAnno_df$start + 1L
centroAnno_df$mid <- floor((centroAnno_df$start + centroAnno_df$end) / 2)
centroAnno_df$display_width <- pmax(centroAnno_df$width, CENTROANNO_DISPLAY_MIN)
centroAnno_df$plot_start <- floor(centroAnno_df$mid - centroAnno_df$display_width / 2)
centroAnno_df$plot_end <- ceiling(centroAnno_df$mid + centroAnno_df$display_width / 2)
centroAnno_df$chr_len <- chr_lengths[centroAnno_df$chrom]
centroAnno_df$plot_start <- pmax(1L, centroAnno_df$plot_start)
centroAnno_df$plot_end <- pmin(centroAnno_df$chr_len, centroAnno_df$plot_end)
centroAnno_gr <- toGRanges(data.frame(
  chr = centroAnno_df$chrom, start = centroAnno_df$plot_start, end = centroAnno_df$plot_end
))

# ----------------------------------------------------------------------------
# Colors
# ----------------------------------------------------------------------------
CENPA_COLOR       <- "#1F5A94"
CENTROANNO_COLOR  <- "#4DAF4A"
CENTROANNO_BORDER <- "#2D6F2A"
TRIANGLE_COLOR    <- "#F2C94C"
IDEOGRAM_FILL     <- "#E5E5E5"
IDEOGRAM_BORDER   <- "#4D4D4D"

# ----------------------------------------------------------------------------
# Plot parameters (plot.type = 4: 1 above, 1 below ideogram)
# ----------------------------------------------------------------------------
pp <- getDefaultPlotParams(plot.type = 4)

pp$leftmargin   <- 0.090
pp$rightmargin  <- 0.025
pp$topmargin    <- 24
pp$bottommargin <- 20

pp$data1height   <- 110   # CENP-A signal (above)
pp$ideogramheight <- 8
pp$data2height   <- 50   # centroAnno intervals (below)

pp$chromosome.margin <- 0.018

# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
FIGURE_WIDTH  <- 8.5
FIGURE_HEIGHT <- 11
FIGURE_RES    <- 300
BASE_PTSIZE   <- 12

output_pdf <- file.path(BASE_DIR, "genome_wide_CENPA_karyotype_all_centroAnno.pdf")
output_png <- file.path(BASE_DIR, "genome_wide_CENPA_karyotype_all_centroAnno_preview.png")

# PNG preview
png(output_png, width = FIGURE_WIDTH, height = FIGURE_HEIGHT,
    units = "in", res = FIGURE_RES, pointsize = BASE_PTSIZE, bg = "white")

par(bg = "white", mar = c(0, 0, 0, 0), oma = c(0, 0, 5, 0), xpd = NA)

kp <- plotKaryotype(
  genome      = custom_genome,
  chromosomes = CHROMOSOMES,
  plot.type   = 4,
  plot.params = pp,
  ideogram.plotter = NULL,
  labels.plotter   = NULL
)

# Ideograms
kpRect(kp, data = custom_genome, y0 = 0, y1 = 1, data.panel = "ideogram",
       col = IDEOGRAM_FILL, border = IDEOGRAM_BORDER, lwd = 0.70)
kpAddChromosomeNames(kp, chr.names = CHROMOSOMES, cex = 0.95)

# --------------------------------------------------------------------------
# Track 1 (above ideogram): CENP-A signal bars
# --------------------------------------------------------------------------
SIGNAL_R0 <- 0.035
SIGNAL_R1 <- 0.60

kpDataBackground(kp, data.panel = 1, r0 = 0.01, r1 = 0.99, color = "white")

for (chrom_i in CHROMOSOMES) {
  chr_gr <- bins_gr[as.character(seqnames(bins_gr)) == chrom_i]
  if (length(chr_gr) == 0) next
  kpBars(kp, data = chr_gr, y0 = 0, y1 = chr_gr$y, data.panel = 1,
         r0 = SIGNAL_R0, r1 = SIGNAL_R1,
         ymin = 0, ymax = chr_max_lookup[chrom_i],
         col = CENPA_COLOR, border = NA)
}

# Triangle markers at highest-mean 100-kb window
for (i in seq_len(nrow(top_regions))) {
  chr_len <- chr_lengths[top_regions$chrom[i]]
  half_width <- min(2200000L, chr_len) / 2
  tip_x <- top_regions$peak_midpoint[i]
  left  <- max(1, tip_x - half_width)
  right <- min(chr_len, tip_x + half_width)
  tip_y <- SIGNAL_R0 + top_regions$peak_fraction[i] * (SIGNAL_R1 - SIGNAL_R0)
  base_y <- tip_y + 0.34

  kpPolygon(kp, chr = top_regions$chrom[i],
            x = c(left, right, tip_x),
            y = c(base_y, base_y, tip_y),
            data.panel = 1, r0 = 0, r1 = 1, ymin = 0, ymax = 1,
            col = TRIANGLE_COLOR, border = NA)
}

# --------------------------------------------------------------------------
# Track 2 (below ideogram): ALL centroAnno satellite arrays
# --------------------------------------------------------------------------
kpDataBackground(kp, data.panel = 2, r0 = 0.02, r1 = 0.98, color = "white")

kpPlotRegions(kp, data = centroAnno_gr, data.panel = 2,
              r0 = 0.08, r1 = 0.94,
              col = CENTROANNO_COLOR, border = CENTROANNO_BORDER, lwd = 0.25)

# Label
kpAddLabels(kp, labels = sprintf("centroAnno\nsatellite arrays\n(n=%d)", nrow(centroAnno_raw)),
            data.panel = 2, r0 = 0.08, r1 = 0.94,
            cex = 0.8, col = CENTROANNO_COLOR, font = 2, pos = 4,
            label.margin = 0.02)

# Title
mtext("Genome-wide CENP-A signal and centroAnno satellite arrays",
      side = 3, outer = TRUE, line = 1.5, adj = 0.5, cex = 1.3, font = 2, col = "black")

# Legend
par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE, xpd = NA)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
legend(x = 0.36, y = 0.10, xjust = 0, yjust = 1,
       legend = c(
         "CENP-A CUT&Tag rep1 (25-kb bins)",
         "Tallest bin in highest-mean 100-kb window",
         sprintf("centroAnno satellite arrays (n=%d)", nrow(centroAnno_raw))
       ),
       pch = c(15, 25, 15),
       col = c(CENPA_COLOR, TRIANGLE_COLOR, CENTROANNO_COLOR),
       pt.bg = c(CENPA_COLOR, TRIANGLE_COLOR, CENTROANNO_COLOR),
       pt.cex = c(1.9, 1.9, 1.9),
       text.font = 1, bty = "n", cex = 1.1,
       x.intersp = 0.85, y.intersp = 1.15)

dev.off()

# PDF
cairo_pdf(output_pdf, width = FIGURE_WIDTH, height = FIGURE_HEIGHT,
          pointsize = BASE_PTSIZE, family = "sans", bg = "white")

par(bg = "white", mar = c(0, 0, 0, 0), oma = c(0, 0, 5, 0), xpd = NA)

kp <- plotKaryotype(
  genome      = custom_genome,
  chromosomes = CHROMOSOMES,
  plot.type   = 4,
  plot.params = pp,
  ideogram.plotter = NULL,
  labels.plotter   = NULL
)

kpRect(kp, data = custom_genome, y0 = 0, y1 = 1, data.panel = "ideogram",
       col = IDEOGRAM_FILL, border = IDEOGRAM_BORDER, lwd = 0.70)
kpAddChromosomeNames(kp, chr.names = CHROMOSOMES, cex = 0.95)

for (chrom_i in CHROMOSOMES) {
  chr_gr <- bins_gr[as.character(seqnames(bins_gr)) == chrom_i]
  if (length(chr_gr) == 0) next
  kpBars(kp, data = chr_gr, y0 = 0, y1 = chr_gr$y, data.panel = 1,
         r0 = SIGNAL_R0, r1 = SIGNAL_R1,
         ymin = 0, ymax = chr_max_lookup[chrom_i],
         col = CENPA_COLOR, border = NA)
}

for (i in seq_len(nrow(top_regions))) {
  chr_len <- chr_lengths[top_regions$chrom[i]]
  half_width <- min(2200000L, chr_len) / 2
  tip_x <- top_regions$peak_midpoint[i]
  left  <- max(1, tip_x - half_width)
  right <- min(chr_len, tip_x + half_width)
  tip_y <- SIGNAL_R0 + top_regions$peak_fraction[i] * (SIGNAL_R1 - SIGNAL_R0)
  base_y <- tip_y + 0.34
  kpPolygon(kp, chr = top_regions$chrom[i],
            x = c(left, right, tip_x), y = c(base_y, base_y, tip_y),
            data.panel = 1, r0 = 0, r1 = 1, ymin = 0, ymax = 1,
            col = TRIANGLE_COLOR, border = NA)
}

kpDataBackground(kp, data.panel = 2, r0 = 0.02, r1 = 0.98, color = "white")
kpPlotRegions(kp, data = centroAnno_gr, data.panel = 2,
              r0 = 0.08, r1 = 0.94,
              col = CENTROANNO_COLOR, border = CENTROANNO_BORDER, lwd = 0.25)

kpAddLabels(kp, labels = sprintf("centroAnno\nsatellite arrays\n(n=%d)", nrow(centroAnno_raw)),
            data.panel = 2, r0 = 0.08, r1 = 0.94,
            cex = 0.8, col = CENTROANNO_COLOR, font = 2, pos = 4,
            label.margin = 0.02)

mtext("Genome-wide CENP-A signal and centroAnno satellite arrays",
      side = 3, outer = TRUE, line = 1.5, adj = 0.5, cex = 1.3, font = 2, col = "black")

par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE, xpd = NA)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
legend(x = 0.36, y = 0.10, xjust = 0, yjust = 1,
       legend = c(
         "CENP-A CUT&Tag rep1 (25-kb bins)",
         "Tallest bin in highest-mean 100-kb window",
         sprintf("centroAnno satellite arrays (n=%d)", nrow(centroAnno_raw))
       ),
       pch = c(15, 25, 15),
       col = c(CENPA_COLOR, TRIANGLE_COLOR, CENTROANNO_COLOR),
       pt.bg = c(CENPA_COLOR, TRIANGLE_COLOR, CENTROANNO_COLOR),
       pt.cex = c(1.9, 1.9, 1.9),
       text.font = 1, bty = "n", cex = 1.1,
       x.intersp = 0.85, y.intersp = 1.15)

dev.off()

cat(sprintf("\nSaved:\n  %s\n  %s\n", output_pdf, output_png))
cat("Done.\n")
