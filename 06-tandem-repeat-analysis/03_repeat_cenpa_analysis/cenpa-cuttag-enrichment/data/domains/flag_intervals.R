#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

# Args
args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]

# Read counts
counts_150 <- fread(file.path(base_dir, "XG_150_setA_counts.txt"),
                     header = FALSE, col.names = c("chrom", "start", "end", "id", "count"))
counts_151 <- fread(file.path(base_dir, "XG_151_setA_counts.txt"),
                     header = FALSE, col.names = c("chrom", "start", "end", "id", "count"))

# Read library sizes
lib_150 <- as.numeric(strsplit(readLines(file.path(base_dir, "XG_150_lib_size.txt")), " ")[[1]][2])
lib_151 <- as.numeric(strsplit(readLines(file.path(base_dir, "XG_151_lib_size.txt")), " ")[[1]][2])

# Compute interval length and normalized CUT&Tag signal (use ilen to avoid collision with base::length)
counts_150[, ilen := end - start]
counts_150[, signal_150 := (count / (ilen / 1000)) / (lib_150 / 1e6)]
counts_151[, ilen := end - start]
counts_151[, signal_151 := (count / (ilen / 1000)) / (lib_151 / 1e6)]

# Merge
dt <- merge(counts_150[, .(chrom, start, end, id, ilen, count_150 = count, signal_150)],
            counts_151[, .(chrom, start, end, ilen, count_151 = count, signal_151)],
            by = c("chrom", "start", "end"))

# Pseudocount
all_signal <- c(dt$signal_150, dt$signal_151)
min_nonzero <- min(all_signal[all_signal > 0], na.rm = TRUE)
pseudocount <- min_nonzero / 2

dt[, log2_signal_150 := log2(signal_150 + pseudocount)]
dt[, log2_signal_151 := log2(signal_151 + pseudocount)]
dt[, mean_cenpa := (log2_signal_150 + log2_signal_151) / 2]

# Per-chromosome median CENP-A
chr_medians <- dt[, .(chr_median = median(mean_cenpa, na.rm = TRUE)), by = chrom]

# Flag intervals above chromosome median
dt <- merge(dt, chr_medians, by = "chrom")
dt[, cenpa_positive := mean_cenpa > chr_median]

# Write flagged intervals BED (cenpa_positive = TRUE)
flagged <- dt[(cenpa_positive), .(chrom, start, end, id)]
fwrite(flagged, file.path(base_dir, "cenpa_positive_intervals.bed"),
       sep = "\t", col.names = FALSE)

# Write full diagnostic table
fwrite(dt, file.path(base_dir, "interval_cenpa_signal.csv"))

# Summary
cat(sprintf("Total intervals: %d\n", nrow(dt)))
cat(sprintf("CENP-A-positive intervals: %d (%.1f%%)\n",
            nrow(flagged), 100 * nrow(flagged) / nrow(dt)))
cat(sprintf("Pseudocount: %.6f\n", pseudocount))
cat(sprintf("Mean CENP-A range: %.3f to %.3f\n", min(dt$mean_cenpa), max(dt$mean_cenpa)))

# Per-chromosome summary
chr_summary <- dt[, .(
    n_total = .N,
    n_positive = sum(cenpa_positive),
    pct_positive = round(100 * sum(cenpa_positive) / .N, 1),
    median_cenpa = round(median(mean_cenpa), 3)
), by = chrom][order(chrom)]
cat("\nPer-chromosome flagging:\n")
print(chr_summary, row.names = FALSE)
