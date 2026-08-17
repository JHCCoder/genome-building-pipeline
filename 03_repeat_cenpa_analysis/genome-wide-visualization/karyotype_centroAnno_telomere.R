# ============================================================================
# Genome-wide CENP-A + centroAnno HOR arrays + TTAGGG telomeric repeats
#
#   Blue bars:
#       CENP-A CUT&Tag rep1 signal (25-kb mean bins, aggregated to 250-kb
#       display bins). Per-chromosome scaling.
#
#   Yellow triangles:
#       Tallest 25-kb bin within the highest-mean 100-kb window.
#
#   Brown ticks:
#       Midpoint locations of centroAnno HOR arrays (from HORs.bed).
#
#   Red ticks (beneath the HOR track):
#       TTAGGG / CCCTAA telomeric repeats (>= 500 bp; exact hexamer scan
#       merged with mismatch-tolerant TRF query).
#
#   Layout:
#       Base-R one-column format. All 30 chromosomes stacked with alternating
#       gray/white row backgrounds. Legend floats over unused right-hand space.
#
#   Self-contained — runs inside local({}) to avoid overwriting globals.
# ============================================================================

local({

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

# ----------------------------------------------------------------------------
# Input paths
# ----------------------------------------------------------------------------

FIGURE_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/genome-wide-repeat-visualization-track"

FAI_FILE <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/degu-genome-browser-pythonVersion/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta.fai"

BW_FILE <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/degu-genome-browser-pythonVersion/XG_150.all.bw"

CENTROANNO_HOR_BED <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-centraAnno/hifiasm-0414/cautils-chrOnly/HORs.bed"

TELOMERE_BED <- file.path(
  FIGURE_DIR,
  "output/telomere_arrays_combined_v2_300bp.bed"
)

ARRAY349_BED <- file.path(
  FIGURE_DIR,
  "output/bin6_349bp_arrays.bed"
)

# ----------------------------------------------------------------------------
# Chromosomes to plot (all 30)
# ----------------------------------------------------------------------------

plot_chromosomes <- c(
  paste0("chr", 1:28), "chrX", "chrY"
)

N_CHROMOSOMES <- length(plot_chromosomes)

# ----------------------------------------------------------------------------
# Read chromosome sizes
# ----------------------------------------------------------------------------

chr_sizes_all <- fread(
  FAI_FILE,
  header = FALSE,
  select = c(1, 2),
  col.names = c("chrom", "size")
)

chr_sizes_all <- chr_sizes_all[chrom %in% plot_chromosomes]
chr_sizes_all[, chromosome_order := match(chrom, plot_chromosomes)]
setorder(chr_sizes_all, chromosome_order)
chr_sizes_all[, chromosome_order := NULL]

chr_lengths <- setNames(
  as.numeric(chr_sizes_all$size),
  chr_sizes_all$chrom
)

max_chromosome_length <- max(chr_lengths, na.rm = TRUE)

# ----------------------------------------------------------------------------
# Load CENP-A CUT&Tag BigWig signal (25-kb mean bins)
# ----------------------------------------------------------------------------

ANALYTICAL_BIN_WIDTH <- 25000L   # 25 kb
DISPLAY_BIN_WIDTH    <- 250000L  # 250 kb (for aggregation)

cat("\nLoading CENP-A BigWig signal ...\n")

bw <- BigWigFile(BW_FILE)

all_bins <- vector("list", N_CHROMOSOMES)
names(all_bins) <- plot_chromosomes

for (chrom_i in plot_chromosomes) {

  chr_len <- chr_lengths[[chrom_i]]

  n_bins <- max(1L, as.integer(ceiling(chr_len / ANALYTICAL_BIN_WIDTH)))

  summarized_signal <- summary(
    bw,
    size  = n_bins,
    type  = "mean",
    which = GRanges(
      seqnames = chrom_i,
      ranges   = IRanges(start = 1, end = chr_len)
    )
  )[[1]]

  all_bins[[chrom_i]] <- data.table(
    chrom = chrom_i,
    start = start(summarized_signal),
    end   = end(summarized_signal),
    score = as.numeric(summarized_signal$score)
  )
}

bins_dt <- rbindlist(all_bins, use.names = TRUE, fill = TRUE)
bins_dt[is.na(score) | !is.finite(score), score := 0]
bins_dt[, plot_score := pmax(score, 0)]

cat(sprintf("Analytical bins loaded: %d\n", nrow(bins_dt)))

# ----------------------------------------------------------------------------
# Aggregate analytical bins to display bins (250 kb)
# ----------------------------------------------------------------------------

bins_dt[
  ,
  display_start := ((start - 1L) %/% DISPLAY_BIN_WIDTH) * DISPLAY_BIN_WIDTH + 1L
]

bins_dt[
  ,
  display_end := pmin(
    display_start + DISPLAY_BIN_WIDTH - 1L,
    as.integer(chr_lengths[chrom])
  )
]

signal_display_dt <- bins_dt[
  ,
  .(signal = mean(plot_score, na.rm = TRUE)),
  by = .(chrom, display_start, display_end)
]

signal_display_dt <- signal_display_dt[is.finite(signal)]

# Per-chromosome max for independent scaling
chromosome_signal_max <- signal_display_dt[
  ,
  .(signal_max = max(signal, na.rm = TRUE)),
  by = chrom
]

chromosome_signal_max[
  !is.finite(signal_max) | signal_max <= 0,
  signal_max := 1
]

chromosome_signal_max <- setNames(
  chromosome_signal_max$signal_max,
  chromosome_signal_max$chrom
)

cat(sprintf(
  "Display bins after aggregation: %d\n",
  nrow(signal_display_dt)
))

# ----------------------------------------------------------------------------
# Find peak CENP-A signal for yellow triangle marker (per chromosome)
# (tallest 25-kb bin within the highest-mean 100-kb window)
# ----------------------------------------------------------------------------

PEAK_WINDOW_BP <- 100000L   # 100-kb sliding window
PEAK_WINDOW_BINS <- as.integer(PEAK_WINDOW_BP / ANALYTICAL_BIN_WIDTH)  # 4 bins

peak_midpoints <- list()

for (chrom_i in plot_chromosomes) {

  chr_bins <- bins_dt[chrom == chrom_i][order(start)]

  if (nrow(chr_bins) >= PEAK_WINDOW_BINS) {
    n_windows <- nrow(chr_bins) - PEAK_WINDOW_BINS + 1L
    window_means <- numeric(n_windows)
    for (w in seq_len(n_windows)) {
      window_means[w] <- mean(
        chr_bins$plot_score[w:(w + PEAK_WINDOW_BINS - 1L)], na.rm = TRUE
      )
    }
    best_window_idx <- which.max(window_means)

    win_start <- best_window_idx
    win_end   <- best_window_idx + PEAK_WINDOW_BINS - 1L
    best_bin_idx <- win_start - 1L + which.max(chr_bins$plot_score[win_start:win_end])

    peak_midpoints[[chrom_i]] <- as.integer(floor(
      (chr_bins$start[best_bin_idx] + chr_bins$end[best_bin_idx]) / 2
    ))
  } else {
    peak_midpoints[[chrom_i]] <- NULL
  }
}

# ----------------------------------------------------------------------------
# Load centroAnno HOR arrays (HORs.bed)
# ----------------------------------------------------------------------------

cat("\nLoading centroAnno HOR arrays ...\n")

hor_raw <- fread(
  CENTROANNO_HOR_BED,
  header = FALSE,
  col.names = c("chrom", "name", "start", "end", "n_monomers", "hor_length", "span_length")
)

hor_raw <- hor_raw[chrom %in% plot_chromosomes]

hor_dt <- hor_raw[, .(
  chrom,
  start,
  end,
  midpoint = as.integer(floor((start + end) / 2))
)]

hor_dt[, midpoint := pmax(1L, pmin(
  as.integer(chr_lengths[chrom]), as.integer(midpoint)
))]

cat(sprintf(
  "centroAnno HOR arrays loaded: %d regions across %d chromosomes\n",
  nrow(hor_dt),
  uniqueN(hor_dt$chrom)
))

# ----------------------------------------------------------------------------
# Load TTAGGG telomeric repeats (midpoints)
# ----------------------------------------------------------------------------

cat("\nLoading TTAGGG telomeric repeats ...\n")

tel_raw <- fread(
  TELOMERE_BED,
  header = FALSE,
  select = c(1, 2, 3),
  col.names = c("chrom", "start", "end")
)

tel_raw <- tel_raw[chrom %in% plot_chromosomes]

tel_dt <- tel_raw[, .(
  chrom,
  start,
  end,
  midpoint = as.integer(floor((start + end) / 2))
)]

tel_dt[, midpoint := pmax(1L, pmin(
  as.integer(chr_lengths[chrom]), as.integer(midpoint)
))]

cat(sprintf(
  "TTAGGG telomeric repeats loaded: %d arrays across %d chromosomes\n",
  nrow(tel_dt),
  uniqueN(tel_dt$chrom)
))

# ----------------------------------------------------------------------------
# Load 349-bp satellite arrays (midpoints)
# ----------------------------------------------------------------------------

cat("\nLoading 349-bp satellite arrays ...\n")

arr349_raw <- fread(
  ARRAY349_BED,
  header = FALSE,
  select = c(1, 2, 3),
  col.names = c("chrom", "start", "end")
)

arr349_raw <- arr349_raw[chrom %in% plot_chromosomes]

arr349_dt <- arr349_raw[, .(
  chrom,
  start,
  end,
  midpoint = as.integer(floor((start + end) / 2))
)]

arr349_dt[, midpoint := pmax(1L, pmin(
  as.integer(chr_lengths[chrom]), as.integer(midpoint)
))]

cat(sprintf(
  "349-bp satellite arrays loaded: %d arrays across %d chromosomes\n",
  nrow(arr349_dt),
  uniqueN(arr349_dt$chrom)
))

# ----------------------------------------------------------------------------
# Figure settings — tall one-column genome-wide layout
# ----------------------------------------------------------------------------

FIGURE_WIDTH_IN  <- 8.5
FIGURE_HEIGHT_IN <- 18
FIGURE_RES_DPI   <- 400
BASE_POINTSIZE   <- 11

# ----------------------------------------------------------------------------
# Colors
# ----------------------------------------------------------------------------

CENPA_COLOR      <- "#2166AC"   # blue bars (CENP-A signal)
HOR_TICK_COLOR   <- "#8B4513"   # brown ticks (centroAnno HOR arrays)
ARR349_TICK_COLOR <- "#1B7837"  # green ticks (349-bp satellite arrays)
TEL_TICK_COLOR   <- "#D43F3A"   # red ticks (TTAGGG telomeres)
TOP_REGION_COLOR <- "#F2C94C"   # yellow peak triangle

IDEOGRAM_FILL    <- "#BDBDBD"
IDEOGRAM_BORDER  <- "#4D4D4D"

# Alternating gray-white chromosome-row layering
ALTERNATING_ROW_FILL <- "#F4F4F4"

# ----------------------------------------------------------------------------
# Text sizes
# ----------------------------------------------------------------------------

TITLE_CEX            <- 1.20
CHROMOSOME_LABEL_CEX <- 0.85
LEGEND_CEX           <- 0.95

# ----------------------------------------------------------------------------
# Physical track dimensions (inches -> NDC)
# ----------------------------------------------------------------------------

SIGNAL_MAX_HEIGHT_IN <- 0.24

IDEOGRAM_HALF_HEIGHT_IN <- 0.014

# Three non-overlapping lower tracks, each ~0.05 in tall (inches below row centre):
#   centroAnno HORs  : 0.030 - 0.080
#   349-bp satellite : 0.080 - 0.130
#   TTAGGG telomeres : 0.130 - 0.180   (kept lowest, no overlap)
HOR_TICK_BOTTOM_OFFSET_IN <- 0.080
HOR_TICK_TOP_OFFSET_IN    <- 0.030
HOR_TICK_LWD              <- 0.70

ARR349_TICK_BOTTOM_OFFSET_IN <- 0.130
ARR349_TICK_TOP_OFFSET_IN    <- 0.080
ARR349_TICK_LWD              <- 0.70

TEL_TICK_BOTTOM_OFFSET_IN <- 0.180
TEL_TICK_TOP_OFFSET_IN    <- 0.130
TEL_TICK_LWD              <- 1.10

SIGNAL_MAX_HEIGHT <- SIGNAL_MAX_HEIGHT_IN / FIGURE_HEIGHT_IN

IDEOGRAM_HALF_HEIGHT <- IDEOGRAM_HALF_HEIGHT_IN / FIGURE_HEIGHT_IN

HOR_TICK_BOTTOM_OFFSET <- HOR_TICK_BOTTOM_OFFSET_IN / FIGURE_HEIGHT_IN
HOR_TICK_TOP_OFFSET    <- HOR_TICK_TOP_OFFSET_IN    / FIGURE_HEIGHT_IN

ARR349_TICK_BOTTOM_OFFSET <- ARR349_TICK_BOTTOM_OFFSET_IN / FIGURE_HEIGHT_IN
ARR349_TICK_TOP_OFFSET    <- ARR349_TICK_TOP_OFFSET_IN    / FIGURE_HEIGHT_IN

TEL_TICK_BOTTOM_OFFSET <- TEL_TICK_BOTTOM_OFFSET_IN / FIGURE_HEIGHT_IN
TEL_TICK_TOP_OFFSET    <- TEL_TICK_TOP_OFFSET_IN    / FIGURE_HEIGHT_IN

# ----------------------------------------------------------------------------
# Triangle settings (inverted yellow peak marker)
# ----------------------------------------------------------------------------

TRIANGLE_HALF_WIDTH_BP <- 1000000L

TRIANGLE_BASE_OFFSET_IN <- 0.270
TRIANGLE_TIP_OFFSET_IN  <- 0.195

TRIANGLE_BASE_OFFSET <- TRIANGLE_BASE_OFFSET_IN / FIGURE_HEIGHT_IN
TRIANGLE_TIP_OFFSET  <- TRIANGLE_TIP_OFFSET_IN  / FIGURE_HEIGHT_IN

# ----------------------------------------------------------------------------
# Fixed layout coordinates
# ----------------------------------------------------------------------------

TITLE_Y <- 0.982

CHROMOSOME_ROW_TOP    <- 0.960
CHROMOSOME_ROW_BOTTOM <- 0.032

CHROMOSOME_LABEL_X <- 0.058
TRACK_START_X      <- 0.080
TRACK_END_X        <- 0.968

# Compute row positions — evenly spaced from top to bottom
row_positions <- seq(
  from = CHROMOSOME_ROW_TOP,
  to   = CHROMOSOME_ROW_BOTTOM,
  length.out = N_CHROMOSOMES
)

# ----------------------------------------------------------------------------
# Floating legend — positioned in unused right-hand space
# ----------------------------------------------------------------------------

# Place legend lower, beside the short chromosomes chr26-chr28 (rows 26-28)
# so it sits in free space to the right of the short ideograms and does not
# overlap chromosome rows.
legend_anchor_row <- row_positions[27]  # beside chr27 (short chromosome)

LEGEND_SYMBOL_X <- 0.360
LEGEND_TEXT_X   <- 0.383

LEGEND_ITEM_GAP_IN <- 0.24
LEGEND_ITEM_GAP <- LEGEND_ITEM_GAP_IN / FIGURE_HEIGHT_IN

LEGEND_CENTER_Y <- legend_anchor_row

LEGEND_ROW_1_Y <- LEGEND_CENTER_Y + 2.0 * LEGEND_ITEM_GAP
LEGEND_ROW_2_Y <- LEGEND_CENTER_Y + 1.0 * LEGEND_ITEM_GAP
LEGEND_ROW_3_Y <- LEGEND_CENTER_Y
LEGEND_ROW_4_Y <- LEGEND_CENTER_Y - 1.0 * LEGEND_ITEM_GAP
LEGEND_ROW_5_Y <- LEGEND_CENTER_Y - 2.0 * LEGEND_ITEM_GAP

# ----------------------------------------------------------------------------
# Title and output filenames
# ----------------------------------------------------------------------------

MAIN_TITLE <- "Genome-wide CENP-A, centroAnno HORs, 349-bp satellite, and TTAGGG telomeres"

OUTPUT_PNG <- file.path(FIGURE_DIR, "genome_wide_CENPA_centroAnno_HORs_telomere.png")
OUTPUT_PDF <- file.path(FIGURE_DIR, "genome_wide_CENPA_centroAnno_HORs_telomere.pdf")
OUTPUT_SVG <- file.path(FIGURE_DIR, "genome_wide_CENPA_centroAnno_HORs_telomere.svg")

# ----------------------------------------------------------------------------
# Coordinate-conversion helper
# ----------------------------------------------------------------------------

genomic_position_to_x <- function(genomic_position) {
  TRACK_START_X +
    (as.numeric(genomic_position) / max_chromosome_length) *
    (TRACK_END_X - TRACK_START_X)
}

# ----------------------------------------------------------------------------
# Draw one chromosome row
# ----------------------------------------------------------------------------

draw_chromosome_row <- function(chrom_i, row_i, row_y) {

  chromosome_length_i <- chr_lengths[chrom_i]
  chromosome_end_x <- genomic_position_to_x(chromosome_length_i)

  # ---- Alternating gray-white row background ----

  if (row_i %% 2L == 1L) {
    row_half_height <- (row_positions[1] - row_positions[2]) / 2 * 0.95
    rect(
      xleft   = 0.005,
      ybottom = row_y - row_half_height,
      xright  = 1.0,
      ytop    = row_y + row_half_height,
      col     = ALTERNATING_ROW_FILL,
      border  = NA
    )
  }

  # ---- Chromosome label ----

  text(
    x       = CHROMOSOME_LABEL_X,
    y       = row_y,
    labels  = chrom_i,
    adj     = c(1, 0.5),
    cex     = CHROMOSOME_LABEL_CEX,
    font    = 1,
    family  = "sans",
    col     = "black"
  )

  # ---- CENP-A signal bars (above ideogram) ----

  chromosome_signal <- signal_display_dt[chrom == chrom_i]

  signal_baseline_y <- row_y + 0.0008

  if (nrow(chromosome_signal) > 0L) {

    signal_max_i <- unname(chromosome_signal_max[chrom_i])

    if (length(signal_max_i) != 1L ||
        !is.finite(signal_max_i) ||
        signal_max_i <= 0) {
      signal_max_i <- 1
    }

    normalized_signal <- pmax(0, pmin(1,
      chromosome_signal$signal / signal_max_i
    ))

    signal_x_left  <- genomic_position_to_x(
      chromosome_signal$display_start - 1L
    )
    signal_x_right <- genomic_position_to_x(
      chromosome_signal$display_end
    )

    signal_y_top <- signal_baseline_y +
      normalized_signal * SIGNAL_MAX_HEIGHT

    rect(
      xleft   = signal_x_left,
      ybottom = signal_baseline_y,
      xright  = signal_x_right,
      ytop    = signal_y_top,
      col     = CENPA_COLOR,
      border  = NA
    )
  }

  # ---- Yellow peak triangle (inverted, tip at peak signal top) ----

  peak_mp <- peak_midpoints[[chrom_i]]

  if (!is.null(peak_mp)) {

    peak_x <- genomic_position_to_x(peak_mp)

    triangle_half_width_x <- (
      TRIANGLE_HALF_WIDTH_BP /
        max_chromosome_length
    ) * (TRACK_END_X - TRACK_START_X)

    triangle_left_x <- max(
      TRACK_START_X,
      peak_x - triangle_half_width_x
    )

    triangle_right_x <- min(
      chromosome_end_x,
      peak_x + triangle_half_width_x
    )

    polygon(
      x = c(triangle_left_x, triangle_right_x, peak_x),
      y = c(
        row_y + TRIANGLE_BASE_OFFSET,
        row_y + TRIANGLE_BASE_OFFSET,
        row_y + TRIANGLE_TIP_OFFSET
      ),
      col     = TOP_REGION_COLOR,
      border  = NA
    )
  }

  # ---- Chromosome ideogram ----

  rect(
    xleft   = TRACK_START_X,
    ybottom = row_y - IDEOGRAM_HALF_HEIGHT,
    xright  = chromosome_end_x,
    ytop    = row_y + IDEOGRAM_HALF_HEIGHT,
    col     = IDEOGRAM_FILL,
    border  = IDEOGRAM_BORDER,
    lwd     = 0.30
  )

  # ---- centroAnno HOR array midpoint ticks ----

  chromosome_hor <- hor_dt[chrom == chrom_i]

  if (nrow(chromosome_hor) > 0L) {

    hor_x <- genomic_position_to_x(chromosome_hor$midpoint)

    segments(
      x0    = hor_x,
      y0    = row_y - HOR_TICK_BOTTOM_OFFSET,
      x1    = hor_x,
      y1    = row_y - HOR_TICK_TOP_OFFSET,
      col   = HOR_TICK_COLOR,
      lwd   = HOR_TICK_LWD,
      lend  = "butt"
    )
  }

  # ---- 349-bp satellite array midpoint ticks (green, beneath the HOR track) ----

  chromosome_arr349 <- arr349_dt[chrom == chrom_i]

  if (nrow(chromosome_arr349) > 0L) {

    arr349_x <- genomic_position_to_x(chromosome_arr349$midpoint)

    segments(
      x0    = arr349_x,
      y0    = row_y - ARR349_TICK_BOTTOM_OFFSET,
      x1    = arr349_x,
      y1    = row_y - ARR349_TICK_TOP_OFFSET,
      col   = ARR349_TICK_COLOR,
      lwd   = ARR349_TICK_LWD,
      lend  = "butt"
    )
  }

  # ---- TTAGGG telomeric repeat ticks (red, lowest track) ----

  chromosome_tel <- tel_dt[chrom == chrom_i]

  if (nrow(chromosome_tel) > 0L) {

    tel_x <- genomic_position_to_x(chromosome_tel$midpoint)

    segments(
      x0    = tel_x,
      y0    = row_y - TEL_TICK_BOTTOM_OFFSET,
      x1    = tel_x,
      y1    = row_y - TEL_TICK_TOP_OFFSET,
      col   = TEL_TICK_COLOR,
      lwd   = TEL_TICK_LWD,
      lend  = "butt"
    )
  }

  invisible(NULL)
}

# ----------------------------------------------------------------------------
# Legend symbols
# ----------------------------------------------------------------------------

draw_square_legend_symbol <- function(x, y, color) {
  symbol_half_width  <- 0.0055
  symbol_half_height <- 0.055 / FIGURE_HEIGHT_IN

  rect(
    xleft   = x - symbol_half_width,
    ybottom = y - symbol_half_height,
    xright  = x + symbol_half_width,
    ytop    = y + symbol_half_height,
    col     = color,
    border  = NA
  )
}

draw_tick_legend_symbol <- function(x, y, color, lwd = 1.7) {
  tick_half_height <- 0.090 / FIGURE_HEIGHT_IN

  segments(
    x0    = x,
    y0    = y - tick_half_height,
    x1    = x,
    y1    = y + tick_half_height,
    col   = color,
    lwd   = lwd,
    lend  = "butt"
  )
}

draw_triangle_legend_symbol <- function(x, y, color) {
  tri_half_width  <- 0.007
  tri_half_height <- 0.070 / FIGURE_HEIGHT_IN

  polygon(
    x = c(x - tri_half_width, x + tri_half_width, x),
    y = c(y + tri_half_height, y + tri_half_height, y - tri_half_height),
    col     = color,
    border  = NA
  )
}

# ----------------------------------------------------------------------------
# Floating legend
# Triangle, CENP-A, HORs, Telomeres
# ----------------------------------------------------------------------------

draw_floating_legend <- function() {

  # Peak CENP-A marker
  draw_triangle_legend_symbol(
    x     = LEGEND_SYMBOL_X,
    y     = LEGEND_ROW_1_Y,
    color = TOP_REGION_COLOR
  )

  text(
    x       = LEGEND_TEXT_X,
    y       = LEGEND_ROW_1_Y,
    labels  = "Peak CENP-A (tallest 25-kb bin in highest-mean 100-kb window)",
    adj     = c(0, 0.5),
    cex     = LEGEND_CEX,
    family  = "sans",
    col     = "black"
  )

  # CENP-A signal
  draw_square_legend_symbol(
    x     = LEGEND_SYMBOL_X,
    y     = LEGEND_ROW_2_Y,
    color = CENPA_COLOR
  )

  text(
    x       = LEGEND_TEXT_X,
    y       = LEGEND_ROW_2_Y,
    labels  = sprintf(
      "CENP-A CUT&Tag rep1 (%s-kb display bins)",
      format(DISPLAY_BIN_WIDTH / 1000, scientific = FALSE, trim = TRUE)
    ),
    adj     = c(0, 0.5),
    cex     = LEGEND_CEX,
    family  = "sans",
    col     = "black"
  )

  # centroAnno HOR arrays
  draw_tick_legend_symbol(
    x     = LEGEND_SYMBOL_X,
    y     = LEGEND_ROW_3_Y,
    color = HOR_TICK_COLOR
  )

  text(
    x       = LEGEND_TEXT_X,
    y       = LEGEND_ROW_3_Y,
    labels  = sprintf(
      "centroAnno HOR arrays (midpoint regions; n = %s)",
      format(nrow(hor_dt), big.mark = ",")
    ),
    adj     = c(0, 0.5),
    cex     = LEGEND_CEX,
    family  = "sans",
    col     = "black"
  )

  # 349-bp satellite arrays (from TRF, period 348-349 bp)
  draw_tick_legend_symbol(
    x     = LEGEND_SYMBOL_X,
    y     = LEGEND_ROW_4_Y,
    color = ARR349_TICK_COLOR
  )

  text(
    x       = LEGEND_TEXT_X,
    y       = LEGEND_ROW_4_Y,
    labels  = sprintf(
      "349-bp satellite arrays (TRF, period 348–349 bp; n = %s)",
      format(nrow(arr349_dt), big.mark = ",")
    ),
    adj     = c(0, 0.5),
    cex     = LEGEND_CEX,
    family  = "sans",
    col     = "black"
  )

  # TTAGGG telomeric repeats
  draw_tick_legend_symbol(
    x     = LEGEND_SYMBOL_X,
    y     = LEGEND_ROW_5_Y,
    color = TEL_TICK_COLOR,
    lwd   = 2.2
  )

  text(
    x       = LEGEND_TEXT_X,
    y       = LEGEND_ROW_5_Y,
    labels  = sprintf(
      "TTAGGG/CCCTAA telomeric repeats ≥ 300 bp (n = %s)",
      format(nrow(tel_dt), big.mark = ",")
    ),
    adj     = c(0, 0.5),
    cex     = LEGEND_CEX,
    family  = "sans",
    col     = "black"
  )

  invisible(NULL)
}

# ----------------------------------------------------------------------------
# Complete figure
# ----------------------------------------------------------------------------

draw_complete_figure <- function() {

  par(
    mar    = c(0, 0, 0, 0),
    oma    = c(0, 0, 0, 0),
    bg     = "white",
    fg     = "black",
    family = "sans",
    xaxs   = "i",
    yaxs   = "i",
    xpd    = NA
  )

  plot.new()

  plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i"
  )

  # Title
  text(
    x       = 0.5,
    y       = TITLE_Y,
    labels  = MAIN_TITLE,
    cex     = TITLE_CEX,
    font    = 2,
    family  = "sans",
    col     = "black",
    adj     = c(0.5, 0.5)
  )

  # Draw all chromosome rows
  for (i in seq_len(N_CHROMOSOMES)) {
    draw_chromosome_row(
      chrom_i = plot_chromosomes[i],
      row_i   = i,
      row_y   = row_positions[i]
    )
  }

  # Floating legend
  draw_floating_legend()

  invisible(NULL)
}

# ----------------------------------------------------------------------------
# Device helpers
# ----------------------------------------------------------------------------

open_png_device <- function(filename) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      filename   = filename,
      width      = FIGURE_WIDTH_IN,
      height     = FIGURE_HEIGHT_IN,
      units      = "in",
      res        = FIGURE_RES_DPI,
      pointsize  = BASE_POINTSIZE,
      background = "white"
    )
  } else if (capabilities("cairo")) {
    png(
      filename   = filename,
      width      = FIGURE_WIDTH_IN,
      height     = FIGURE_HEIGHT_IN,
      units      = "in",
      res        = FIGURE_RES_DPI,
      pointsize  = BASE_POINTSIZE,
      type       = "cairo",
      bg         = "white"
    )
  } else {
    png(
      filename   = filename,
      width      = FIGURE_WIDTH_IN,
      height     = FIGURE_HEIGHT_IN,
      units      = "in",
      res        = FIGURE_RES_DPI,
      pointsize  = BASE_POINTSIZE,
      bg         = "white"
    )
  }
}

open_pdf_device <- function(filename) {
  if (capabilities("cairo")) {
    cairo_pdf(
      filename   = filename,
      width      = FIGURE_WIDTH_IN,
      height     = FIGURE_HEIGHT_IN,
      pointsize  = BASE_POINTSIZE,
      family     = "sans",
      bg         = "white"
    )
  } else {
    pdf(
      file        = filename,
      width       = FIGURE_WIDTH_IN,
      height      = FIGURE_HEIGHT_IN,
      pointsize   = BASE_POINTSIZE,
      family      = "sans",
      useDingbats = FALSE,
      bg          = "white"
    )
  }
}

open_svg_device <- function(filename) {
  svg(
    filename   = filename,
    width      = FIGURE_WIDTH_IN,
    height     = FIGURE_HEIGHT_IN,
    pointsize  = BASE_POINTSIZE,
    bg         = "white"
  )
}

# ----------------------------------------------------------------------------
# Render
# ----------------------------------------------------------------------------

cat("\nDrawing genome-wide CENP-A + centroAnno HORs + TTAGGG telomeres ...\n")

open_png_device(OUTPUT_PNG)
draw_complete_figure()
invisible(dev.off())

cat(sprintf("Saved PNG to: %s\n", normalizePath(OUTPUT_PNG, mustWork = FALSE)))

open_pdf_device(OUTPUT_PDF)
draw_complete_figure()
invisible(dev.off())

cat(sprintf("Saved PDF to: %s\n", normalizePath(OUTPUT_PDF, mustWork = FALSE)))

open_svg_device(OUTPUT_SVG)
draw_complete_figure()
invisible(dev.off())

cat(sprintf("Saved SVG to: %s\n", normalizePath(OUTPUT_SVG, mustWork = FALSE)))

cat("\nDone.\n")

})  # End of local()
