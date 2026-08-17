# ============================================================================
# CENP-A signal + HiCAT HOR arrays — chr4 and chr25
#
#   Blue bars:
#       CENP-A CUT&Tag rep1 signal (25-kb mean bins, aggregated to 250-kb
#       display bins). Per-chromosome scaling.
#
#   Green ticks:
#       Midpoint locations of all HiCAT HOR arrays (>= 100 monomers each).
#       All monomer types shown as a single track — no type separation.
#
#   Layout:
#       Base-R one-column format adapted from the centroAnno genome-wide
#       overview. Stacked chr4 and chr25 rows with alternating gray/white
#       row backgrounds.
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
# Chromosomes to plot (top to bottom)
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

cat(sprintf("chr4:  %s bp\n", format(chr_lengths["chr4"],  big.mark = ",")))
cat(sprintf("chr25: %s bp\n", format(chr_lengths["chr25"], big.mark = ",")))

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
  "Display bins after aggregation: %d (chr4: %d, chr25: %d)\n",
  nrow(signal_display_dt),
  signal_display_dt[chrom == "chr4", .N],
  signal_display_dt[chrom == "chr25", .N]
))

# ----------------------------------------------------------------------------
# Load HiCAT HOR arrays (all monomer types merged, one color only)
# ----------------------------------------------------------------------------

hor_raw_dt <- rbindlist(
  lapply(plot_chromosomes, function(chrom) {
    fpath <- file.path(HICAT_BED_DIR, paste0("hicat_", chrom, "_hor_arrays.bed"))
    dt <- fread(
      fpath,
      skip = 3,
      col.names = c(
        "chrom", "start", "end",
        "n_monomers", "dominant_type", "n_types"
      )
    )
    return(dt)
  }),
  use.names = TRUE,
  fill      = TRUE
)

# Keep only substantial arrays (>= 100 monomers)
hor_sub_dt <- hor_raw_dt[n_monomers >= 100]

# Merge nearby arrays (all types) into contiguous HOR-dense regions
HOR_MERGE_GAP <- 200000L
setorder(hor_sub_dt, chrom, start)

merged_list <- list()
for (chrom_i in plot_chromosomes) {
  sub <- hor_sub_dt[chrom == chrom_i]
  if (nrow(sub) == 0) next

  cur_start <- sub$start[1]
  cur_end   <- sub$end[1]

  if (nrow(sub) > 1) {
    for (r in 2:nrow(sub)) {
      if (sub$start[r] - cur_end <= HOR_MERGE_GAP) {
        cur_end <- max(cur_end, sub$end[r])
      } else {
        merged_list[[length(merged_list) + 1]] <- list(
          chrom = chrom_i, start = cur_start, end = cur_end
        )
        cur_start <- sub$start[r]
        cur_end   <- sub$end[r]
      }
    }
  }
  merged_list[[length(merged_list) + 1]] <- list(
    chrom = chrom_i, start = cur_start, end = cur_end
  )
}

hor_regions_dt <- rbindlist(merged_list)

# Compute midpoints for tick marks
hor_regions_dt[, midpoint := floor((start + end) / 2)]
hor_regions_dt[, midpoint := pmax(1L, pmin(
  as.integer(chr_lengths[chrom]), as.integer(midpoint)
))]

cat(sprintf(
  "HOR arrays: %d raw -> %d merged regions (all types combined)\n",
  nrow(hor_sub_dt),
  nrow(hor_regions_dt)
))

# ----------------------------------------------------------------------------
# Figure settings
# ----------------------------------------------------------------------------

FIGURE_WIDTH_IN  <- 8.5
FIGURE_HEIGHT_IN <- 2.6
FIGURE_RES_DPI   <- 400
BASE_POINTSIZE   <- 11

# ----------------------------------------------------------------------------
# Colors — matching the centroAnno genome-wide overview
# ----------------------------------------------------------------------------

CENPA_COLOR      <- "#2166AC"   # blue bars (CENP-A signal)
HOR_TICK_COLOR   <- "#238B45"   # green ticks (HiCAT HOR arrays)

IDEOGRAM_FILL    <- "#BDBDBD"
IDEOGRAM_BORDER  <- "#4D4D4D"

ALTERNATING_ROW_FILL <- "#F4F4F4"

# ----------------------------------------------------------------------------
# Text sizes
# ----------------------------------------------------------------------------

TITLE_CEX            <- 1.10
CHROMOSOME_LABEL_CEX <- 0.90
LEGEND_CEX           <- 0.80

# ----------------------------------------------------------------------------
# Physical track dimensions (inches -> NDC)
# ----------------------------------------------------------------------------

SIGNAL_MAX_HEIGHT_IN <- 0.22

IDEOGRAM_HALF_HEIGHT_IN <- 0.020

HOR_TICK_BOTTOM_OFFSET_IN <- 0.165
HOR_TICK_TOP_OFFSET_IN    <- 0.055
HOR_TICK_LWD              <- 0.75

SIGNAL_MAX_HEIGHT <- SIGNAL_MAX_HEIGHT_IN / FIGURE_HEIGHT_IN

IDEOGRAM_HALF_HEIGHT <- IDEOGRAM_HALF_HEIGHT_IN / FIGURE_HEIGHT_IN

HOR_TICK_BOTTOM_OFFSET <- HOR_TICK_BOTTOM_OFFSET_IN / FIGURE_HEIGHT_IN
HOR_TICK_TOP_OFFSET    <- HOR_TICK_TOP_OFFSET_IN    / FIGURE_HEIGHT_IN

# ----------------------------------------------------------------------------
# Floating legend
# ----------------------------------------------------------------------------

LEGEND_SYMBOL_X <- 0.390
LEGEND_TEXT_X   <- 0.411

LEGEND_ITEM_GAP_IN <- 0.22
LEGEND_ITEM_GAP <- LEGEND_ITEM_GAP_IN / FIGURE_HEIGHT_IN

LEGEND_CENTER_Y <- 0.60

LEGEND_ROW_1_Y <- LEGEND_CENTER_Y + LEGEND_ITEM_GAP
LEGEND_ROW_2_Y <- LEGEND_CENTER_Y

# ----------------------------------------------------------------------------
# Title and output filenames
# ----------------------------------------------------------------------------

MAIN_TITLE <- "CENP-A signal and HiCAT HOR arrays -- chr4 and chr25"

OUTPUT_PNG <- file.path(
  HICAT_BED_DIR,
  "hicat_chr4_chr25_one_column_preview.png"
)

OUTPUT_PDF <- file.path(
  HICAT_BED_DIR,
  "hicat_chr4_chr25_one_column.pdf"
)

# ----------------------------------------------------------------------------
# Fixed layout coordinates
# ----------------------------------------------------------------------------

TITLE_Y <- 0.960

CHROMOSOME_ROW_TOP    <- 0.900
CHROMOSOME_ROW_BOTTOM <- 0.080

CHROMOSOME_LABEL_X <- 0.070
TRACK_START_X      <- 0.105
TRACK_END_X        <- 0.965

# ----------------------------------------------------------------------------
# Chromosome-row Y positions
# ----------------------------------------------------------------------------

chromosome_row_y <- setNames(
  seq(
    from       = CHROMOSOME_ROW_TOP,
    to         = CHROMOSOME_ROW_BOTTOM,
    length.out = length(plot_chromosomes)
  ),
  plot_chromosomes
)

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

  # ---- Alternating gray/white row background ----

  if (row_i %% 2L == 1L) {
    rect(
      xleft   = CHROMOSOME_LABEL_X - 0.050,
      ybottom = row_y - 0.018,
      xright  = TRACK_END_X,
      ytop    = row_y + 0.026,
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

  signal_baseline_y <- row_y + 0.0035

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

draw_tick_legend_symbol <- function(x, y, color) {
  tick_half_height <- 0.090 / FIGURE_HEIGHT_IN

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

# ----------------------------------------------------------------------------
# Floating legend
# ----------------------------------------------------------------------------

draw_floating_legend <- function() {

  # CENP-A signal
  draw_square_legend_symbol(
    x     = LEGEND_SYMBOL_X,
    y     = LEGEND_ROW_1_Y,
    color = CENPA_COLOR
  )

  text(
    x       = LEGEND_TEXT_X,
    y       = LEGEND_ROW_1_Y,
    labels  = sprintf(
      "CENP-A CUT&Tag rep1 (%s-kb mean display bins)",
      format(DISPLAY_BIN_WIDTH / 1000, scientific = FALSE, trim = TRUE)
    ),
    adj     = c(0, 0.5),
    cex     = LEGEND_CEX,
    family  = "sans",
    col     = "black"
  )

  # HiCAT HOR arrays (all types combined)
  draw_tick_legend_symbol(
    x     = LEGEND_SYMBOL_X,
    y     = LEGEND_ROW_2_Y,
    color = HOR_TICK_COLOR
  )

  text(
    x       = LEGEND_TEXT_X,
    y       = LEGEND_ROW_2_Y,
    labels  = sprintf(
      "HiCAT HOR arrays (midpoint ticks; n = %s regions)",
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

  # Draw chromosome rows
  for (row_i in seq_along(plot_chromosomes)) {
    chrom_i <- plot_chromosomes[row_i]
    draw_chromosome_row(
      chrom_i = chrom_i,
      row_i   = row_i,
      row_y   = unname(chromosome_row_y[chrom_i])
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

# ----------------------------------------------------------------------------
# Render
# ----------------------------------------------------------------------------

cat("\nDrawing CENP-A + HiCAT overview ...\n")

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

cat("\nDone.\n")

})  # End of local()
