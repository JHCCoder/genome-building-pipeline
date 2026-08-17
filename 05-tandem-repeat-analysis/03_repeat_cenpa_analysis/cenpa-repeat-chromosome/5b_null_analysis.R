#!/usr/bin/env Rscript
# ============================================================================
# 5b_null_analysis.R -- compare observed band composition vs matched-shuffle null
#
# For each repeat family, the "concentration" metric is the fraction of the
# family's bp that overlaps each distance band from the CENP-A domain
# (core / core->250kb / 250kb->1Mb / remainder). Bands are constructed per
# chromosome from the CENP-A domain (same geometry as Step 4). Overlap is
# computed by actual bp intersection (not midpoint), so large arrays spanning
# multiple bands are handled correctly.
#
# Observed: fraction of real-array bp per band, one value per chromosome.
# Null:     fraction of matched-shuffle placement bp per band (same geometry).
# We report, per family x band: observed_mean, null_mean, null_2.5/97.5, z, p.
#
# Usage: Rscript 5b_null_analysis.R <work_dir>
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
DATA <- file.path(WORK, "data")
RES  <- file.path(WORK, "results")
SRC  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
dir.create(RES, showWarnings = FALSE)

BAND_NEAR <- 250000
BAND_FAR  <- 1000000

sizes <- fread(file.path(SRC, "data", "chrom_sizes.txt"),
               header = FALSE, col.names = c("chrom", "len"))
sizes <- sizes[chrom %in% paste0("chr", c(1:28, "X", "Y"))]
domains <- fread(file.path(DATA, "domains", "cenpa_domains_weighted.csv"))
setkey(domains, chrom)

# ---- band construction (same as Step 4) ----
make_bands <- function(core_start, core_end, chrom_len) {
  core <- data.table(start = core_start, end = core_end, band = "core")
  near <- rbind(
    data.table(start = max(0, core_start - BAND_NEAR), end = core_start, band = "core_to_250kb"),
    data.table(start = core_end, end = min(chrom_len, core_end + BAND_NEAR), band = "core_to_250kb"))
  far <- rbind(
    data.table(start = max(0, core_start - BAND_FAR), end = max(0, core_start - BAND_NEAR), band = "250kb_to_1Mb"),
    data.table(start = min(chrom_len, core_end + BAND_NEAR), end = min(chrom_len, core_end + BAND_FAR), band = "250kb_to_1Mb"))
  masked <- rbind(core, near, far)[order(start)]
  # remainder = chromosome minus union of core/near/far
  rem <- list(); cur <- 0
  for (i in seq_len(nrow(masked))) {
    s <- max(cur, masked[i, start]); e <- masked[i, end]
    if (s > cur) rem[[length(rem) + 1]] <- data.table(start = cur, end = s, band = "remainder")
    cur <- max(cur, e)
  }
  if (cur < chrom_len) rem[[length(rem) + 1]] <- data.table(start = cur, end = chrom_len, band = "remainder")
  rbindlist(list(core, near, far, rbindlist(rem)))
}

# band overlap for a set of intervals [s,e) -> data.table(chrom,start,end,...)
band_frac <- function(iv, chr) {
  if (!(chr %in% domains$chrom)) return(NULL)
  dom <- domains[domains$chrom == chr]
  clen <- sizes$len[sizes$chrom == chr]
  bands <- make_bands(dom$core_start, dom$core_end, clen)
  iv <- iv[iv$chrom == chr]
  tot <- sum(iv[, end - start])
  out <- list()
  for (i in seq_len(nrow(bands))) {
    bs <- bands$start[i]; be <- bands$end[i]; bn <- bands$band[i]
    ivb <- iv[iv$end > bs & iv$start < be, ]
    ov <- if (nrow(ivb) == 0) 0 else sum(pmax(0, pmin(ivb$end, be) - pmax(ivb$start, bs)))
    out[[i]] <- data.table(family = iv$family[1], chrom = chr, band = bn,
                           frac = if (tot > 0) ov / tot else NA_real_)
  }
  rbindlist(out)
}

# ---- observed ----
real <- list()
fam_name <- c(bin6_349 = "349 bp", bin4_195 = "195 bp", bin8_389 = "389 bp")
for (fam in names(fam_name)) {
  bin <- sub("_.*", "", fam)
  bed <- fread(file.path(SRC, "period-enrichment", "data", "merged", "arrays",
                         paste0(bin, "_arrays.bed")),
               header = FALSE, col.names = c("chrom", "start", "end", "id", "cnt"))
  bed[, family := fam]
  for (chr in unique(bed$chrom)) {
    fb <- band_frac(bed, chr)
    if (!is.null(fb)) real[[length(real) + 1]] <- fb
  }
}
real_dt <- rbindlist(real)
real_dt[is.na(frac), frac := 0]
real_dt[, family := fam_name[family]]

# ---- null placements ----
null <- fread(file.path(RES, "matched_null_placements.csv"))
null[, family := fcase(family == "bin6_349", "349 bp",
                       family == "bin4_195", "195 bp",
                       family == "bin8_389", "389 bp")]
null_dt <- list()
for (chr in unique(null$chrom)) {
  if (!(chr %in% domains$chrom)) next
  dom <- domains[domains$chrom == chr]
  clen <- sizes$len[sizes$chrom == chr]
  bands <- make_bands(dom$core_start, dom$core_end, clen)
  for (fam_n in unique(null$family)) {
    nc <- null[null$chrom == chr & null$family == fam_n, ]
    if (nrow(nc) == 0) next
    tot_all <- sum(nc$null_end - nc$null_start)
    for (i in seq_len(nrow(bands))) {
      bs <- bands$start[i]; be <- bands$end[i]; bn <- bands$band[i]
      ovb <- nc[nc$null_end > bs & nc$null_start < be, ]
      ov <- if (nrow(ovb) == 0) 0 else sum(pmax(0, pmin(ovb$null_end, be) - pmax(ovb$null_start, bs)))
      null_dt[[length(null_dt) + 1]] <- data.table(
        family = fam_n, chrom = chr, band = bn,
        frac = if (tot_all > 0) ov / tot_all else NA_real_)
    }
  }
}
null_dt <- rbindlist(null_dt)
null_dt[is.na(frac), frac := 0]

# ---- compare ----
comp <- real_dt[, .(obs_mean = mean(frac), obs_sd = sd(frac), n_chr = .N),
                by = .(family, band)]
null_agg <- null_dt[, .(null_mean = mean(frac), null_sd = sd(frac)), by = .(family, band)]
out <- merge(comp, null_agg, by = c("family", "band"))
out[, z := (obs_mean - null_mean) / ifelse(null_sd == 0, NA, null_sd)]
out[, p_greater := pnorm(z, lower.tail = FALSE)]
out[, p_two := 2 * pmin(p_greater, 1 - p_greater)]
setorder(out, family, band)
fwrite(out, file.path(RES, "matched_null_band_enrichment.csv"))
cat("=== observed vs matched-null, per family x band (fraction of family bp) ===\n")
print(out, row.names = FALSE)
cat("DONE 5b_null_analysis.R\n")
