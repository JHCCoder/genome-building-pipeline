#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]

# Read merged domains
domains <- fread(file.path(base_dir, "domains_with_interval_counts.bed"),
                 header = FALSE, col.names = c("chrom", "start", "end", "domain_id", "size", "n_intervals"))

# Read CENP-A signal (for mean per domain)
signal <- fread(file.path(base_dir, "interval_cenpa_signal.csv"))

# Count flagged intervals per original centroAnno interval
flagged <- fread(file.path(base_dir, "cenpa_positive_intervals.bed"),
                 header = FALSE, col.names = c("chrom", "start", "end", "id"))

cat("=== Domain QC Summary ===\n\n")
cat(sprintf("Total merged domains: %d\n", nrow(domains)))
cat(sprintf("Median domain size: %.1f kb\n", median(domains$size) / 1000))
cat(sprintf("Mean domain size: %.1f kb\n", mean(domains$size) / 1000))
cat(sprintf("Max domain size: %.1f kb\n", max(domains$size) / 1000))
cat(sprintf("Singletons (1 interval): %d (%.1f%%)\n",
            sum(domains$n_intervals == 1),
            100 * sum(domains$n_intervals == 1) / nrow(domains)))
cat(sprintf("Total flagged intervals in domains: %d / %d\n",
            sum(domains$n_intervals), nrow(flagged)))

cat("\nPer-chromosome domain counts:\n")
chr_counts <- domains[, .(n_domains = .N,
                           total_size_mb = round(sum(size) / 1e6, 2),
                           median_size_kb = round(median(size) / 1000, 1),
                           total_intervals = sum(n_intervals)), by = chrom]
chr_counts <- chr_counts[order(chrom)]
print(chr_counts, row.names = FALSE)

# Size distribution
cat("\nDomain size distribution:\n")
size_breaks <- c(0, 10e3, 50e3, 100e3, 250e3, 500e3, 1e6, 5e6, Inf)
size_labels <- c("<10kb", "10-50kb", "50-100kb", "100-250kb", "250-500kb", "500kb-1Mb", "1-5Mb", ">5Mb")
domains[, size_bin := cut(size, breaks = size_breaks, labels = size_labels, right = FALSE)]
size_dist <- domains[, .N, by = size_bin][order(size_bin)]
print(size_dist, row.names = FALSE)

# chr1 check
chr1_domains <- domains[chrom == "chr1"]
cat(sprintf("\nchr1: %d domains (was 1 super-domain with pure distance merge at 500kb)\n", nrow(chr1_domains)))
cat(sprintf("  Largest chr1 domain: %.1f kb, %d intervals\n",
            max(chr1_domains$size) / 1000, chr1_domains[which.max(size), n_intervals]))

# Fraction of flagged intervals per chromosome
cat("\nFraction of centroAnno intervals flagged as CENP-A-positive:\n")
chr_flag <- signal[, .(n_total = .N, n_flagged = sum(cenpa_positive)), by = chrom]
chr_flag[, pct := round(100 * n_flagged / n_total, 1)]
chr_flag <- chr_flag[order(chrom)]
print(chr_flag, row.names = FALSE)

# Write QC table
fwrite(chr_counts, file.path(base_dir, "domains_per_chromosome.csv"))
fwrite(size_dist, file.path(base_dir, "domain_size_distribution.csv"))
