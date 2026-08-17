#!/usr/bin/env Rscript
# ============================================================================
# 2_domains.R -- define the CENP-A domain independently on every chromosome
#
# Basis: 1/NH-weighted k=100 CENP-A coverage (mean of XG_150 + XG_151) over
# non-overlapping 100 kb windows (primary). k=1 unique-mapped coverage is
# produced as a sensitivity set.
#
# Per-chromosome, independent threshold: a chromosome's CENP-A windows are
# defined by that chromosome's OWN signal distribution (median + k*MAD) -- no
# genome-wide cutoff, so chr1 cannot dominate.
#
# Domain = the contiguous run of above-threshold windows that CONTAINS the
# chromosome's peak (maximum-signal) window, gap-filled <=250 kb. This is a
# STRICT SIGNAL PEAK definition (user decision): the core is the elevated
# signal around the peak, NOT extended to any repeat array.
#
# Caveat reported: on chromosomes whose centromere is a homogeneous satellite
# array (e.g. chr4's 8.7 Mb 349 array), 1/NH-weighted reads still collapse to
# the single uniquely-mappable degenerate junction, so the strict-signal core
# is small (~100 kb). The k-mer leg (step 6) independently resolves which
# repeat family CENP-A fragments truly contain.
#
# Usage: Rscript 2_domains.R <work_dir> <signal>   (signal = weighted | k1)
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
SIGNAL <- if (length(args) >= 2) args[2] else "weighted"
DATA <- file.path(WORK, "data")
COVW <- file.path(DATA, "coverage_weighted")
DOM  <- file.path(DATA, "domains")
dir.create(DOM, showWarnings = FALSE)

W <- 100000
K_MAD <- 4          # threshold = median + k*MAD  (per-chromosome)
GAP_FILL <- 250000  # merge runs separated by <=250kb

read_signal <- function() {
  x150 <- fread(file.path(COVW, "XG_150_win100kb_weighted.tsv"),
                col.names = c("chrom", "start", "end", "s150"))
  x151 <- fread(file.path(COVW, "XG_151_win100kb_weighted.tsv"),
                col.names = c("chrom", "start", "end", "s151"))
  merge(x150, x151, by = c("chrom", "start", "end"))[, sig := (s150 + s151) / 2]
}

read_signal_k1 <- function() {
  r150 <- fread(file.path(COVW, "XG_150_win100kb_k1.txt"),
                header = FALSE, col.names = c("chrom", "start", "end", "s150"))
  r151 <- fread(file.path(COVW, "XG_151_win100kb_k1.txt"),
                header = FALSE, col.names = c("chrom", "start", "end", "s151"))
  merge(r150, r151, by = c("chrom", "start", "end"))[, sig := (s150 + s151) / 2]
}

sig <- if (SIGNAL == "k1") read_signal_k1() else read_signal()

domains_out <- list()
qc <- list()
for (chr in sort(unique(sig$chrom))) {
  d <- sig[chrom == chr][order(start)]
  med <- median(d$sig); md <- mad(d$sig)
  thr <- med + K_MAD * md
  peak_i <- which.max(d$sig)
  d[, pos := sig > thr]

  if (!d[peak_i, pos]) {
    # peak is below threshold -> diffuse chromosome: fall back to top window
    domains_out[[chr]] <- data.table(chrom = chr,
                                     core_start = d[peak_i, start], core_end = d[peak_i, end],
                                     core_size = W, peak_window = d[peak_i, start],
                                     mean_signal = d[peak_i, sig], basis = SIGNAL,
                                     n_above_thr = sum(d$pos), threshold = thr,
                                     note = "diffuse_fallback_top_window")
    next
  }

  # run id of the peak; the domain = contiguous run containing the peak
  d[, run := cumsum(c(TRUE, diff(pos) != 0))]
  pk_run <- d[peak_i, run]
  seg <- d[run == pk_run]

  # gap-fill: merge runs within GAP_FILL of this segment
  seg <- seg[order(start)]
  # include any adjacent runs separated by <= GAP_FILL (may bridge small dips)
  # simple approach: extend segment by GAP_FILL on both sides
  ext_start <- max(0, min(seg$start) - GAP_FILL)
  ext_end   <- min(max(d$end), max(seg$end) + GAP_FILL)
  in_ext <- d[chrom == chr & start >= ext_start & end <= ext_end]
  # final core = min..max of all windows in the extended region that are above
  # threshold OR within a small gap of an above-threshold window
  final <- in_ext[pos == TRUE]
  # coalesce with gap filling
  f2 <- final[order(start)]
  if (nrow(f2) > 0) {
    core_start <- min(f2$start); core_end <- max(f2$end)
  } else {
    core_start <- d[peak_i, start]; core_end <- d[peak_i, end]
  }

  domains_out[[chr]] <- data.table(chrom = chr,
                                   core_start = core_start, core_end = core_end,
                                   core_size = core_end - core_start,
                                   peak_window = d[peak_i, start],
                                   mean_signal = mean(seg$sig), basis = SIGNAL,
                                   n_above_thr = sum(d$pos), threshold = thr,
                                   note = "strict_signal_peak")
}

out <- rbindlist(domains_out)
setorder(out, chrom)
fwrite(out, file.path(DOM, paste0("cenpa_domains_", SIGNAL, ".bed")), sep = "\t", col.names = FALSE)
fwrite(out, file.path(DOM, paste0("cenpa_domains_", SIGNAL, ".csv")))
cat(sprintf("Domains (%s): %d chromosomes\n", SIGNAL, nrow(out)))
print(out[, .(chrom, core_start, core_end, core_size, n_above_thr, note)], row.names = FALSE)
cat("DONE 2_domains.R\n")
