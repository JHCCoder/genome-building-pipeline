#!/usr/bin/env Rscript
# ============================================================================
# 4_composition_9bin.R -- 9-bin TRF composition by distance from CENP-A domain
#
# For each chromosome: bands = core ; core->250kb ; 250kb->1Mb ; remainder.
# Repeat composition = fraction of each band's sequence covered by each of the
# 9 TRF period bins (bin1..bin9) plus "other" (no TRF call).
#
# Two outputs:
#   1. Per-chromosome (N=30 per bin x band) -- chromosomes as biological units
#   2. Genomewide pooled -- combine ALL chromosomes into one CENP-A peak:
#      sum all core bands, all core->250kb bands, etc., then recompute the
#      9-bin stacked composition over the pooled bp.
#
# Usage: Rscript 4_composition_9bin.R <work_dir>
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
DATA <- file.path(WORK, "data")
RES  <- file.path(WORK, "results")
dir.create(RES, showWarnings = FALSE)

BAND_NEAR <- 250000
BAND_FAR  <- 1000000
BIN_LABELS <- c("1-10 bp","11-50 bp","51-192 bp","193-195 bp","196-347 bp",
                "348-349 bp","350-385 bp","386-390 bp","391+ bp")

# ---- inputs -----------------------------------------------------------------
domains <- fread(file.path(DATA, "domains", "cenpa_domains_weighted.csv"))
sizes <- fread(file.path(DATA, "..", "..", "cenpa-cuttag-enrichment", "data", "chrom_sizes.txt"),
               header = FALSE, col.names = c("chrom", "len"))
sizes <- sizes[chrom %in% domains$chrom]

# per-bin merged TRF beds
per_bin <- list()
for (b in 1:9) {
  f <- file.path(DATA, "repdensity", "per_bin", sprintf("bin%d_merged.bed", b))
  if (!file.exists(f)) stop("missing per-bin bed: ", f)
  d <- fread(f, header = FALSE, col.names = c("chrom", "start", "end"))
  d[, bin_id := b]
  per_bin[[b]] <- d
}
rep9 <- rbindlist(per_bin)
setkey(rep9, chrom, start, end)

# ---- band construction (per chromosome) ------------------------------------
# bands: core ; [core_start-near, core_start) + (core_end, core_end+near] ;
#        [core_start-far, core_start-near) + (core_end+near, core_end+far] ;
#        remainder = chromosome minus the above.
make_bands <- function(core_start, core_end, chrom_len) {
  bands <- list()
  # core (clip to chromosome)
  bands[["core"]] <- data.table(start = core_start, end = core_end)
  # near ring
  bands[["core_to_250kb"]] <- rbind(
    data.table(start = max(0, core_start - BAND_NEAR), end = core_start),
    data.table(start = core_end, end = min(chrom_len, core_end + BAND_NEAR))
  )
  # far ring
  bands[["250kb_to_1Mb"]] <- rbind(
    data.table(start = max(0, core_start - BAND_FAR), end = max(0, core_start - BAND_NEAR)),
    data.table(start = min(chrom_len, core_end + BAND_NEAR), end = min(chrom_len, core_end + BAND_FAR))
  )
  # remainder = whole chromosome minus all three
  masked <- rbind(bands[["core"]], bands[["core_to_250kb"]], bands[["250kb_to_1Mb"]])
  masked <- masked[order(start)]
  # subtract union of masked from [0, chrom_len]
  rem <- data.table(start = 0, end = chrom_len)
  out <- list()
  cur <- 0
  for (i in seq_len(nrow(masked))) {
    s <- max(cur, masked[i, start]); e <- masked[i, end]
    if (s > cur) out[[length(out) + 1]] <- data.table(start = cur, end = s)
    cur <- max(cur, e)
  }
  if (cur < chrom_len) out[[length(out) + 1]] <- data.table(start = cur, end = chrom_len)
  bands[["remainder"]] <- rbindlist(out)
  bands
}

# ---- overlap bp between a bin set and a band --------------------------------
overlap_bp <- function(bin_dt, chr, band) {
  # band is a data.table of [start,end) segments on this chromosome
  total <- 0
  for (i in seq_len(nrow(band))) {
    s <- band[i, start]; e <- band[i, end]
    hits <- bin_dt[bin_dt$chrom == chr & end > s & start < e, ]
    if (nrow(hits) == 0) next
    os <- pmax(hits$start, s); oe <- pmin(hits$end, e)
    total <- total + sum(pmax(0, oe - os))
  }
  total
}

# ---- per-chromosome 9-bin composition ---------------------------------------
rows <- list()
for (chr in domains$chrom) {
  dom <- domains[chrom == chr]
  cs <- dom$core_start; ce <- dom$core_end
  clen <- sizes[chrom == chr, len]
  bands <- make_bands(cs, ce, clen)
  band_bp <- vapply(bands, function(b) sum(b[, end - start]), numeric(1))

  for (bn in names(bands)) {
    band <- bands[[bn]]
    bb <- band_bp[[bn]]
    # per-bin overlap
    for (b in 1:9) {
      ov <- overlap_bp(per_bin[[b]], chr, band)
      rows[[length(rows) + 1]] <- data.table(
        chrom = chr, band = bn, band_bp = bb,
        bin_id = b, bin_label = BIN_LABELS[b],
        repeat_bp = ov, frac_band = if (bb > 0) ov / bb else NA_real_
      )
    }
    # "other" = band bp not covered by any TRF (union of all 9 bins)
    allov <- sum(vapply(1:9, function(b) overlap_bp(per_bin[[b]], chr, band), numeric(1)))
    rows[[length(rows) + 1]] <- data.table(
      chrom = chr, band = bn, band_bp = bb,
      bin_id = 0L, bin_label = "other (no TRF)",
      repeat_bp = pmax(0, bb - allov), frac_band = if (bb > 0) pmax(0, bb - allov) / bb else NA_real_
    )
  }
}
comp <- rbindlist(rows)
comp[is.na(frac_band), frac_band := 0]
fwrite(comp, file.path(RES, "repeat_composition_by_band_9bin.csv"))
cat("Per-chromosome 9-bin composition:", nrow(comp), "rows\n")

# ---- genomewide pooled: combine all chromosomes into one CENP-A peak -------
# Pool each band class across chromosomes, then recompute 9-bin fractions
# over pooled bp (bp-weighted mean of per-chromosome fractions).
gw <- comp[, .(pooled_bp = sum(band_bp), pooled_repeat_bp = sum(repeat_bp)), by = .(band, bin_id, bin_label)]
gw[, frac_band := ifelse(pooled_bp > 0, pooled_repeat_bp / pooled_bp, NA_real_)]
gw[is.na(frac_band), frac_band := 0]
gw[, view := "genomewide_pooled"]
setorder(gw, band, bin_id)
fwrite(gw, file.path(RES, "repeat_composition_genomewide_9bin.csv"))
cat("Genomewide pooled 9-bin composition:", nrow(gw), "rows\n")

# ---- also emit a wide table for plotting -----------------------------------
comp_wide <- dcast(comp[bin_id > 0], chrom + band ~ bin_id, value.var = "frac_band")
setnames(comp_wide, as.character(1:9), paste0("bin", 1:9))
fwrite(comp_wide, file.path(RES, "repeat_composition_by_band_9bin_wide.csv"))
cat("DONE 4_composition_9bin.R\n")
