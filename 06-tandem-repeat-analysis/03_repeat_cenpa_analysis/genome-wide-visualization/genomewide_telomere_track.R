# ============================================================================
# Genome-wide track figure: CENP-A + centroAnno + 349-bp satellite + TTAGGG telomeres
#
# 30 chromosomes (chr1-28, chrX, chrY), one row each.
# Tracks (from top to bottom in each row):
#   blue bars   — CENP-A CUT&Tag rep1 signal (250-kb display bins, per-chrom scaled)
#   green ticks — centroAnno-predicted centromeric repeat regions
#   orange ticks — 349-bp satellite arrays (bin6, bedtools merge -d 0)
#   red ticks   — TTAGGG / CCCTAA telomeric arrays (>= 500 bp, exact + TRF combined)
# Self-contained.
# ============================================================================
local({

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

WDIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/genome-wide-repeat-visualization-track"
FIGDIR <- WDIR

FAI_FILE <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/degu-genome-browser-pythonVersion/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta.fai"
BW_FILE  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/degu-genome-browser-pythonVersion/XG_150.all.bw"

CENTRO_BED  <- file.path(WDIR, "output/centroAnno_repeat_regions.bed")
ARRAY349_BED <- file.path(WDIR, "output/bin6_349bp_arrays.bed")
TELOMERE_BED <- file.path(WDIR, "output/telomere_arrays_combined_500bp.bed")

plot_chromosomes <- c(paste0("chr", 1:28), "chrX", "chrY")
N_CHROMOSOMES <- length(plot_chromosomes)

chr_sizes <- fread(FAI_FILE, header = FALSE, select = c(1, 2), col.names = c("chrom", "size"))
chr_sizes <- chr_sizes[chrom %in% plot_chromosomes]
chr_sizes[, order_ := match(chrom, plot_chromosomes)]
setorder(chr_sizes, order_)
chr_lengths <- setNames(as.numeric(chr_sizes$size), chr_sizes$chrom)
max_chr_len <- max(chr_lengths)

# ---------------------------------------------------------------------------
# CENP-A signal (25-kb analytical bins -> 250-kb display bins)
# ---------------------------------------------------------------------------
ANALYTICAL_BIN <- 25000L
DISPLAY_BIN    <- 250000L
bw <- BigWigFile(BW_FILE)
all_bins <- vector("list", N_CHROMOSOMES)
names(all_bins) <- plot_chromosomes
for (ch in plot_chromosomes) {
  n_bins <- max(1L, as.integer(ceiling(chr_lengths[[ch]] / ANALYTICAL_BIN)))
  s <- summary(bw, size = n_bins, type = "mean",
               which = GRanges(ch, IRanges(1, chr_lengths[[ch]])))[[1]]
  all_bins[[ch]] <- data.table(chrom = ch, start = start(s), end = end(s),
                               score = as.numeric(s$score))
}
bins <- rbindlist(all_bins)
bins[is.na(score) | !is.finite(score), score := 0]
bins[, plot_score := pmax(score, 0)]
bins[, disp_start := ((start - 1L) %/% DISPLAY_BIN) * DISPLAY_BIN + 1L]
bins[, disp_end := pmin(disp_start + DISPLAY_BIN - 1L, as.integer(chr_lengths[chrom]))]
sig <- bins[, .(signal = mean(plot_score, na.rm = TRUE)), by = .(chrom, disp_start, disp_end)]
sig <- sig[is.finite(signal)]
chr_max <- setNames(sig[, .(m = max(signal, na.rm = TRUE)), by = chrom]$m,
                    sig[, .(m = max(signal, na.rm = TRUE)), by = chrom]$chrom)
chr_max[!is.finite(chr_max) | chr_max <= 0] <- 1

# ---------------------------------------------------------------------------
# BED tracks (midpoints for ticks)
# ---------------------------------------------------------------------------
read_mid <- function(path) {
  dt <- fread(path, header = FALSE, select = c(1, 2, 3), col.names = c("chrom", "start", "end"))
  dt <- dt[chrom %in% plot_chromosomes]
  dt[, mid := pmax(1L, pmin(as.integer(chr_lengths[chrom]), as.integer(floor((start + 1L + end) / 2))))]
  dt
}
centro_dt   <- read_mid(CENTRO_BED)
arr349_dt   <- read_mid(ARRAY349_BED)
tel_dt      <- read_mid(TELOMERE_BED)

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
FIG_W <- 8.5
FIG_H <- 18
RES   <- 400
TITLE_Y   <- 0.982
ROW_TOP   <- 0.960
ROW_BOT   <- 0.030
LABEL_X   <- 0.058
TRACK_X0  <- 0.080
TRACK_X1  <- 0.968
row_y <- seq(ROW_TOP, ROW_BOT, length.out = N_CHROMOSOMES)

xpos <- function(bp) TRACK_X0 + (bp / max_chr_len) * (TRACK_X1 - TRACK_X0)

IDEO_HALF <- 0.014 / FIG_H
CENPA_COL <- "#2166AC"
CENTRO_COL <- "#238B45"
ARR349_COL <- "#E6852C"
TEL_COL    <- "#D43F3A"

draw_row <- function(ch, i, y) {
  clen <- chr_lengths[[ch]]
  # alternating row bg
  if (i %% 2L == 1L) {
    rh <- (row_y[1] - row_y[2]) / 2 * 0.95
    rect(0.005, y - rh, 1.0, y + rh, col = "#F4F4F4", border = NA)
  }
  text(LABEL_X, y, labels = ch, adj = c(1, 0.5), cex = 0.85, col = "black")

  # CENP-A signal bars
  csig <- sig[chrom == ch]
  if (nrow(csig) > 0) {
    m <- unname(chr_max[ch]); if (!is.finite(m) || m <= 0) m <- 1
    ns <- pmax(0, pmin(1, csig$signal / m))
    rect(xpos(csig$disp_start - 1L), y + 0.001,
         xpos(csig$disp_end), y + 0.001 + ns * (0.24 / FIG_H),
         col = CENPA_COL, border = NA)
  }
  # ideogram
  rect(TRACK_X0, y - IDEO_HALF, xpos(clen), y + IDEO_HALF,
       col = "#BDBDBD", border = "#4D4D4D", lwd = 0.3)
  # centroAnno ticks (above ideogram)
  cg <- centro_dt[chrom == ch]
  if (nrow(cg) > 0) segments(xpos(cg$mid), y + 0.035, xpos(cg$mid), y + 0.11, col = CENTRO_COL, lwd = 0.7)
  # 349bp ticks (below ideogram, orange)
  a349 <- arr349_dt[chrom == ch]
  if (nrow(a349) > 0) segments(xpos(a349$mid), y - 0.11, xpos(a349$mid), y - 0.035, col = ARR349_COL, lwd = 0.7)
  # telomere ticks (bottom, red, thick)
  tt <- tel_dt[chrom == ch]
  if (nrow(tt) > 0) segments(xpos(tt$mid), y - 0.135, xpos(tt$mid), y - 0.075, col = TEL_COL, lwd = 1.6)
  invisible(NULL)
}

# legend symbols
leg_sq <- function(x, y, col) rect(x - 0.004, y - 0.055/FIG_H, x + 0.004, y + 0.055/FIG_H, col = col, border = NA)
leg_tick <- function(x, y, col, lwd=1.5) segments(x, y - 0.09/FIG_H, x, y + 0.09/FIG_H, col = col, lwd = lwd)

# ---------------------------------------------------------------------------
# render
# ---------------------------------------------------------------------------
render <- function(device = "png", path) {
  if (device == "png") {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(path, width = FIG_W, height = FIG_H, units = "in", res = RES, background = "white")
    } else {
      png(path, width = FIG_W, height = FIG_H, units = "in", res = RES, type = "cairo", bg = "white")
    }
  } else if (device == "pdf") {
    cairo_pdf(path, width = FIG_W, height = FIG_H, family = "sans", bg = "white")
  } else {
    svg(path, width = FIG_W, height = FIG_H, bg = "white")
  }
  par(mar = c(0,0,0,0), oma = c(0,0,0,0), bg = "white", xaxs = "i", yaxs = "i", xpd = NA)
  plot.new(); plot.window(c(0,1), c(0,1), xaxs = "i", yaxs = "i")
  text(0.5, TITLE_Y, "Genome-wide tracks: CENP-A · centroAnno · 349-bp satellite · TTAGGG telomeres",
       cex = 1.15, font = 2, adj = c(0.5, 0.5))
  for (i in seq_len(N_CHROMOSOMES)) draw_row(plot_chromosomes[i], i, row_y[i])

  # legend (floating, right side near bottom rows)
  ly <- row_y[27]
  leg_sq(0.36, ly, CENPA_COL); text(0.385, ly, "CENP-A CUT&Tag rep1 (250-kb bins)", adj = c(0,0.5), cex = 0.9)
  leg_tick(0.36, ly - 0.032, CENTRO_COL); text(0.385, ly - 0.032, sprintf("centroAnno regions (n = %s)", format(nrow(centro_dt), big.mark=",")), adj = c(0,0.5), cex = 0.9)
  leg_tick(0.36, ly - 0.064, ARR349_COL); text(0.385, ly - 0.064, sprintf("349-bp satellite arrays (n = %s)", format(nrow(arr349_dt), big.mark=",")), adj = c(0,0.5), cex = 0.9)
  leg_tick(0.36, ly - 0.096, TEL_COL, lwd=2); text(0.385, ly - 0.096, sprintf("TTAGGG/CCCTAA telomere arrays >=500 bp (n = %s)", format(nrow(tel_dt), big.mark=",")), adj = c(0,0.5), cex = 0.9)
  dev.off()
}

cat(sprintf("n centroAnno=%d, n 349bp=%d, n telomere=%d\n", nrow(centro_dt), nrow(arr349_dt), nrow(tel_dt)))
render("png", file.path(FIGDIR, "genomewide_CENPA_centroAnno_349bp_telomere.png"))
render("pdf", file.path(FIGDIR, "genomewide_CENPA_centroAnno_349bp_telomere.pdf"))
cat("saved.\n")

}) # local
