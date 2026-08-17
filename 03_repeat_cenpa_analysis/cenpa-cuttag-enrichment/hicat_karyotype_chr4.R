# ============================================================================
# CENP-A signal + HiCAT HOR annotations — chr4 + chr25
#
#   Blue bars:
#       CENP-A CUT&Tag rep1 signal (25-kb mean bins, aggregated to 250-kb
#       display bins). Per-chromosome scaling.
#
#   Yellow triangle:
#       Tallest 25-kb bin within the highest-mean 100-kb window.
#
#   Red ticks:
#       Midpoint locations of individual HiCAT HOR calls (repeat_n >= 5).
#       All monomer types shown as a single track — no type separation.
#
#   Layout:
#       Base-R two-row format. Chr4 (top) + chr25 (bottom), each with
#       expanded vertical space for the signal track.
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

HICAT_BED_DIR <- paste0(
  "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/",
  "figure/cenpa-cuttag-enrichment"
)

FAI_FILE <- paste0(
  "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/",
  "degu-genome-browser-pythonVersion/",
  "assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta.fai"
)

BW_FILE <- paste0(
  "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/",
  "degu-genome-browser-pythonVersion/",
  "XG_150.all.bw"
)

# ----------------------------------------------------------------------------
# Chromosomes to plot — chr4 + chr25
# ----------------------------------------------------------------------------

plot_chromosomes <- c("chr4", "chr25")

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

for (chrom_i in plot_chromosomes) {
  cat(sprintf("%s:  %s bp\n", chrom_i, format(chr_lengths[chrom_i], big.mark = ",")))
}

max_chromosome_length <- max(chr_lengths, na.rm = TRUE)

# ----------------------------------------------------------------------------
# Load CENP-A CUT&Tag BigWig signal (25-kb mean bins)
# ----------------------------------------------------------------------------

ANALYTICAL_BIN_WIDTH <- 25000L   # 25 kb
DISPLAY_BIN_WIDTH    <- 250000L  # 250 kb (for aggregation)

cat("\nLoading CENP-A BigWig signal ...\n")

bw <- BigWigFile(BW_FILE)

all_bins <- vector("list", length(plot_chromosomes))
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
    # Slide 100-kb window and compute mean score
    n_windows <- nrow(chr_bins) - PEAK_WINDOW_BINS + 1L
    window_means <- numeric(n_windows)
    for (w in seq_len(n_windows)) {
      window_means[w] <- mean(chr_bins$plot_score[w:(w + PEAK_WINDOW_BINS - 1L)], na.rm = TRUE)
    }
    best_window_idx <- which.max(window_means)

    # Within best window, find the tallest individual 25-kb bin
    win_start <- best_window_idx
    win_end   <- best_window_idx + PEAK_WINDOW_BINS - 1L
    best_bin_idx <- win_start - 1L + which.max(chr_bins$plot_score[win_start:win_end])

    peak_midpoints[[chrom_i]] <- as.integer(floor(
      (chr_bins$start[best_bin_idx] + chr_bins$end[best_bin_idx]) / 2
    ))
    peak_score <- chr_bins$plot_score[best_bin_idx]

    cat(sprintf(
      "Peak CENP-A (%s): bin %d (%d-%d), midpoint = %s, score = %.3f\n",
      chrom_i,
      best_bin_idx,
      chr_bins$start[best_bin_idx],
      chr_bins$end[best_bin_idx],
      format(peak_midpoints[[chrom_i]], big.mark = ","),
      peak_score
    ))
  } else {
    peak_midpoints[[chrom_i]] <- NULL
    cat(sprintf("Not enough bins for peak detection on %s.\n", chrom_i))
  }
}

# ----------------------------------------------------------------------------
# Load HiCAT HOR calls directly from merged out_all_layer.xls
# (individual HOR units, not pre-merged arrays)
# ----------------------------------------------------------------------------

HICAT_BASE_DIR <- paste0(
  "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/",
  "code/command-line-script/genome-annotation/HiCAT/HiCAT_genome"
)

hor_list <- vector("list", length(plot_chromosomes))

for (i in seq_along(plot_chromosomes)) {

  chrom_i <- plot_chromosomes[i]

  merged_dir <- file.path(HICAT_BASE_DIR, paste0(chrom_i, "_merged"))

  hor_raw_dt <- fread(
    file.path(merged_dir, "out_all_layer.xls"),
    header = FALSE,
    col.names = c("hor_id", "start", "end", "repeat_n", "layer", "type")
  )

  hor_raw_dt[, chrom := chrom_i]

  # Keep only HORs with repeat_n >= 10 to filter low-confidence calls
  HOR_MIN_REPEAT <- 10L
  hor_sub_dt <- hor_raw_dt[repeat_n >= 10]

  # Use individual HOR midpoints as tick positions — no merging
  hor_regions_dt <- hor_sub_dt[, .(
    chrom,
    start,
    end,
    midpoint = as.integer(floor((start + end) / 2))
  )]

  hor_regions_dt[, midpoint := pmax(1L, pmin(
    as.integer(chr_lengths[chrom]), as.integer(midpoint)
  ))]

  cat(sprintf(
    "HOR calls (%s): %d total, %d with repeat_n >= %d shown as ticks\n",
    chrom_i,
    nrow(hor_raw_dt),
    nrow(hor_regions_dt),
    HOR_MIN_REPEAT
  ))

  hor_list[[i]] <- hor_regions_dt
}

hor_regions_dt <- rbindlist(hor_list, use.names = TRUE, fill = TRUE)

# ----------------------------------------------------------------------------
# Figure settings — expanded for two-chromosome view
# ----------------------------------------------------------------------------

FIGURE_WIDTH_IN  <- 8.5
FIGURE_HEIGHT_IN <- 3.8
FIGURE_RES_DPI   <- 400
BASE_POINTSIZE   <- 11

# ----------------------------------------------------------------------------
# Colors — matching the centroAnno genome-wide overview
# ----------------------------------------------------------------------------

CENPA_COLOR         <- "#2166AC"   # blue bars (CENP-A signal)
HOR_TICK_COLOR      <- "#DC3220"   # red ticks (HiCAT HOR arrays)
TOP_REGION_COLOR    <- "#F2C94C"   # yellow peak triangle

IDEOGRAM_FILL    <- "#BDBDBD"
IDEOGRAM_BORDER  <- "#4D4D4D"

# Alternating gray-white chromosome-row layering
ALTERNATING_ROW_FILL <- "#F4F4F4"

# ----------------------------------------------------------------------------
# Text sizes
# ----------------------------------------------------------------------------

TITLE_CEX            <- 1.05
CHROMOSOME_LABEL_CEX <- 0.90
LEGEND_CEX           <- 0.75

# ----------------------------------------------------------------------------
# Physical track dimensions (inches -> NDC)
# ----------------------------------------------------------------------------

SIGNAL_MAX_HEIGHT_IN <- 0.30

IDEOGRAM_HALF_HEIGHT_IN <- 0.022

HOR_TICK_BOTTOM_OFFSET_IN <- 0.200
HOR_TICK_TOP_OFFSET_IN    <- 0.065
HOR_TICK_LWD              <- 0.80

SIGNAL_MAX_HEIGHT <- SIGNAL_MAX_HEIGHT_IN / FIGURE_HEIGHT_IN

IDEOGRAM_HALF_HEIGHT <- IDEOGRAM_HALF_HEIGHT_IN / FIGURE_HEIGHT_IN

HOR_TICK_BOTTOM_OFFSET <- HOR_TICK_BOTTOM_OFFSET_IN / FIGURE_HEIGHT_IN
HOR_TICK_TOP_OFFSET    <- HOR_TICK_TOP_OFFSET_IN    / FIGURE_HEIGHT_IN

# ----------------------------------------------------------------------------
# Triangle settings (inverted yellow peak marker)
# ----------------------------------------------------------------------------

TRIANGLE_HALF_WIDTH_BP <- 1000000L

TRIANGLE_BASE_OFFSET_IN <- 0.385
TRIANGLE_TIP_OFFSET_IN  <- 0.315

TRIANGLE_BASE_OFFSET <- TRIANGLE_BASE_OFFSET_IN / FIGURE_HEIGHT_IN
TRIANGLE_TIP_OFFSET  <- TRIANGLE_TIP_OFFSET_IN  / FIGURE_HEIGHT_IN

# ----------------------------------------------------------------------------
# Floating legend
# ----------------------------------------------------------------------------

LEGEND_SYMBOL_X <- 0.530
LEGEND_TEXT_X   <- 0.551

LEGEND_ITEM_GAP_IN <- 0.24
LEGEND_ITEM_GAP <- LEGEND_ITEM_GAP_IN / FIGURE_HEIGHT_IN

LEGEND_CENTER_Y <- 0.185   # below chr25, clear of coordinate labels

LEGEND_ROW_1_Y <- LEGEND_CENTER_Y + 2 * LEGEND_ITEM_GAP
LEGEND_ROW_2_Y <- LEGEND_CENTER_Y + LEGEND_ITEM_GAP
LEGEND_ROW_3_Y <- LEGEND_CENTER_Y

# ----------------------------------------------------------------------------
# Title and output filenames
# ----------------------------------------------------------------------------

MAIN_TITLE <- "Selected chromosomes CENP-A signal and HiCAT HOR annotations"

OUTPUT_PNG <- file.path(
  HICAT_BED_DIR,
  "hicat_chr4_chr25_one_column.png"
)

OUTPUT_PDF <- file.path(
  HICAT_BED_DIR,
  "hicat_chr4_chr25_one_column.pdf"
)

OUTPUT_SVG <- file.path(
  HICAT_BED_DIR,
  "hicat_chr4_chr25_one_column.svg"
)

# ----------------------------------------------------------------------------
# Shared coordinate labels — between chr4 and chr25 rows
# Guidelines extend up into chr4 and down into chr25
# ----------------------------------------------------------------------------

draw_shared_coordinate_labels <- function(row_y_top, row_y_bottom) {

  label_step_mb <- 20  # every 20 Mb
  label_step_bp <- label_step_mb * 1e6

  # Labels up to the longest chromosome
  n_labels <- floor(max_chromosome_length / label_step_bp)

  # Y-positions: label sits midway between the two rows
  gap_label_y  <- (row_y_top + row_y_bottom) / 2

  for (i in seq_len(n_labels)) {
    pos_bp <- i * label_step_bp
    pos_x  <- genomic_position_to_x(pos_bp)

    # Guideline UP to chr4 ideogram (always)
    segments(
      x0  = pos_x,
      y0  = gap_label_y + 0.015,
      x1  = pos_x,
      y1  = row_y_top - IDEOGRAM_HALF_HEIGHT,
      col = "grey80",
      lwd = 0.4,
      lty = "dashed",
      lend = "butt"
    )

    # Guideline DOWN to chr25 ideogram (only if within chr25 bounds)
    if (pos_bp <= chr_lengths["chr25"]) {
      segments(
        x0  = pos_x,
        y0  = gap_label_y - 0.015,
        x1  = pos_x,
        y1  = row_y_bottom + IDEOGRAM_HALF_HEIGHT,
        col = "grey80",
        lwd = 0.4,
        lty = "dashed",
        lend = "butt"
      )
    }

    # Tick mark
    segments(
      x0  = pos_x,
      y0  = gap_label_y - 0.015,
      x1  = pos_x,
      y1  = gap_label_y + 0.015,
      col = "grey60",
      lwd = 0.5,
      lend = "butt"
    )

    # Label
    text(
      x       = pos_x,
      y       = gap_label_y - 0.030,
      labels  = paste0(i * label_step_mb, " Mb"),
      cex     = 0.55,
      font    = 1,
      family  = "sans",
      col     = "grey50",
      adj     = c(0.5, 1)
    )
  }
}

# ----------------------------------------------------------------------------
# Fixed layout coordinates
# ----------------------------------------------------------------------------

TITLE_Y <- 0.930

# Two-row layout: chr4 on top, chr25 on bottom
CHROMOSOME_ROW_TOP_Y    <- 0.700   # chr4
CHROMOSOME_ROW_BOTTOM_Y <- 0.280   # chr25

CHROMOSOME_LABEL_X <- 0.070
TRACK_START_X      <- 0.105
TRACK_END_X        <- 0.965

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

draw_chromosome_row <- function(chrom_i, row_y, row_i) {

  chromosome_length_i <- chr_lengths[chrom_i]
  chromosome_end_x <- genomic_position_to_x(chromosome_length_i)

  # ---- Alternating gray-white row background ----

  if (row_i %% 2L == 1L) {
    rect(
      xleft   = CHROMOSOME_LABEL_X - 0.050,
      ybottom = row_y - 0.195,
      xright  = TRACK_END_X,
      ytop    = row_y + 0.115,
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

  signal_baseline_y <- row_y + 0.004

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
    lwd     = 0.65
  )

  # ---- HiCAT HOR array midpoint ticks (all types, one color) ----

  chromosome_hor <- hor_regions_dt[chrom == chrom_i]

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

  invisible(NULL)
}

# ----------------------------------------------------------------------------
# Legend symbols
# ----------------------------------------------------------------------------

draw_square_legend_symbol <- function(x, y, color) {
  symbol_half_width  <- 0.0060
  symbol_half_height <- 0.070 / FIGURE_HEIGHT_IN

  rect(
    xleft   = x - symbol_half_width,
    ybottom = y - symbol_half_height,
    xright  = x + symbol_half_width,
    ytop    = y + symbol_half_height,
    col     = color,
    border  = NA
  )
}

draw_tick_legend_symbol <- function(x, y, color) {
  tick_half_height <- 0.110 / FIGURE_HEIGHT_IN

  segments(
    x0    = x,
    y0    = y - tick_half_height,
    x1    = x,
    y1    = y + tick_half_height,
    col   = color,
    lwd   = 1.7,
    lend  = "butt"
  )
}

draw_triangle_legend_symbol <- function(x, y, color) {
  tri_half_width  <- 0.008
  tri_half_height <- 0.090 / FIGURE_HEIGHT_IN

  polygon(
    x = c(x - tri_half_width, x + tri_half_width, x),
    y = c(y + tri_half_height, y + tri_half_height, y - tri_half_height),
    col     = color,
    border  = NA
  )
}

# ----------------------------------------------------------------------------
# Floating legend
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

  # HiCAT HOR arrays
  draw_tick_legend_symbol(
    x     = LEGEND_SYMBOL_X,
    y     = LEGEND_ROW_3_Y,
    color = HOR_TICK_COLOR
  )

  text(
    x       = LEGEND_TEXT_X,
    y       = LEGEND_ROW_3_Y,
    labels  = sprintf(
      "HiCAT HORs (repeat_n >= %d; n = %s regions)",
      HOR_MIN_REPEAT,
      format(nrow(hor_regions_dt), big.mark = ",")
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

  # Draw chr4 (top row, row 1 = gray background)
  draw_chromosome_row(chrom_i = "chr4",  row_y = CHROMOSOME_ROW_TOP_Y,    row_i = 1L)
  # Draw chr25 (bottom row, row 2 = white background)
  draw_chromosome_row(chrom_i = "chr25", row_y = CHROMOSOME_ROW_BOTTOM_Y, row_i = 2L)

  # Shared coordinate labels between the two rows
  draw_shared_coordinate_labels(
    row_y_top    = CHROMOSOME_ROW_TOP_Y,
    row_y_bottom = CHROMOSOME_ROW_BOTTOM_Y
  )

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

cat("\nDrawing CENP-A + HiCAT chr4 + chr25 overview ...\n")

open_png_device(OUTPUT_PNG)
draw_complete_figure()
invisible(dev.off())

if (requireNamespace("IRdisplay", quietly = TRUE)) {
  IRdisplay::display_png(file = OUTPUT_PNG)
}

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
