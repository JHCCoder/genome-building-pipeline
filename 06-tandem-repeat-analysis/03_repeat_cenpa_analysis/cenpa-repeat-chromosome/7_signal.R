#!/usr/bin/env Rscript
# ============================================================================
# 7_signal.R -- CENP-A signal per family x band, unique vs 1/NH-weighted
#
# Combines the k=1 and 1/NH-weighted k=100 window coverage (Step 1) with the
# distance bands (Step 4 geometry) to report, per family x band:
#   - mean k=1 CENP-A signal   (unique mappings)
#   - mean 1/NH-weighted signal (multimappers weighted fractionally)
#   - both normalized by band bp and by library size (per 1e6 fragments)
# This is the "unique vs weighted" quantitative comparison in the three-way set.
#
# Usage: Rscript 7_signal.R <work_dir>
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
DATA <- file.path(WORK, "data")
COVW <- file.path(DATA, "coverage_weighted")
RES  <- file.path(WORK, "results")
SRC  <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
dir.create(RES, showWarnings = FALSE)

BAND_NEAR <- 250000
BAND_FAR  <- 1000000

sizes <- fread(file.path(SRC, "data", "chrom_sizes.txt"),
               header = FALSE, col.names = c("chrom", "len"))
sizes <- sizes[chrom %in% paste0("chr", c(1:28, "X", "Y"))]
domains <- fread(file.path(DATA, "domains", "cenpa_domains_weighted.csv"))

# window coverage, mean over CENPA reps (XG_150/151)
win150w <- fread(file.path(COVW, "XG_150_win100kb_weighted.tsv"),
                 col.names = c("chrom", "start", "end", "w150"))
win151w <- fread(file.path(COVW, "XG_151_win100kb_weighted.tsv"),
                 col.names = c("chrom", "start", "end", "w151"))
win150k <- fread(file.path(COVW, "XG_150_win100kb_k1.txt"),
                 header = FALSE, col.names = c("chrom", "start", "end", "k150"))
win151k <- fread(file.path(COVW, "XG_151_win100kb_k1.txt"),
                 header = FALSE, col.names = c("chrom", "start", "end", "k151"))

w <- merge(win150w, win151w, by = c("chrom", "start", "end"))
k <- merge(win150k, win151k, by = c("chrom", "start", "end"))
w[, wmean := (w150 + w151) / 2]
k[, kmean := (k150 + k151) / 2]
win <- merge(w[, .(chrom, start, end, wmean)], k[, .(chrom, start, end, kmean)],
             by = c("chrom", "start", "end"))
win[, mid := (start + end) / 2]

# band for each window center (per chromosome, vectorized)
band_of <- function(pos, cs, ce) {
  in_core <- pos >= cs & pos < ce
  d <- pmin(abs(pos - cs), abs(pos - ce))
  band <- ifelse(in_core, "core",
          ifelse(d <= BAND_NEAR, "core_to_250kb",
          ifelse(d <= BAND_FAR, "250kb_to_1Mb", "remainder")))
  band
}
win[, band := "unassigned"]
for (chr in domains$chrom) {
  win[win$chrom == chr, band := band_of(mid, domains[domains$chrom == chr, core_start],
                                              domains[domains$chrom == chr, core_end])]
}

# normalize signal to per-1e6-fragments-per-bp scale (comparable across reps)
# lib sizes (fragment counts)
lib_size <- c(XG_150 = 27156788, XG_151 = 0)  # filled below
lib_size[["XG_151"]] <- as.numeric(system("wc -l < /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/data/fragments/XG_151_fragments.bed", intern = TRUE))
win[, wcpm := wmean / mean(lib_size) * 1e6]
win[, kcpm := kmean / mean(lib_size) * 1e6]

agg <- win[, .(weighted_mean = mean(wcpm), k1_mean = mean(kcpm),
               n_windows = .N, band_bp = sum(end - start)),
           by = .(chrom, band)]
fwrite(agg, file.path(RES, "cenpa_signal_by_band_window.csv"))

# also per family arrays (349/195/389) x band: mean signal over overlapping windows
fam_out <- list()
for (bin in c("bin6", "bin4", "bin8")) {
  fam <- fread(file.path(SRC, "period-enrichment", "data", "merged", "arrays",
                         paste0(bin, "_arrays.bed")),
               header = FALSE, col.names = c("chrom", "start", "end", "id", "cnt"))
  fam[, family := bin]
  # windows overlapping each array
  setkey(win, chrom, start, end)
  setkey(fam, chrom, start, end)
  ov <- foverlaps(fam, win, nomatch = 0)
  if (nrow(ov) > 0) {
    fam_agg <- ov[, .(n_arrays = uniqueN(id), n_overlap_windows = .N,
                      weighted_mean = mean(wcpm), k1_mean = mean(kcpm)),
                  by = .(family, band)]
  } else {
    fam_agg <- data.table(family = bin, band = character())
  }
  fam_out[[bin]] <- fam_agg
}
fam_res <- rbindlist(fam_out, fill = TRUE)
fwrite(fam_res, file.path(RES, "cenpa_signal_by_band_family.csv"))
cat("=== CENP-A signal by band (windows) ===\n")
print(agg[order(chrom, band)], row.names = FALSE)
cat("DONE 7_signal.R\n")
