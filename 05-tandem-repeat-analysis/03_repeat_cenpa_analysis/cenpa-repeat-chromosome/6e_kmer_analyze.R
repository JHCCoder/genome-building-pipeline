#!/usr/bin/env Rscript
# ============================================================================
# 6e_kmer_analyze.R -- assembly-independent 31-mer content of CENP-A vs control
#
# Question: does CENP-A fragment DNA truly contain 349 / 195 / 389 sequence,
# regardless of where the fragments map? We count family-specific 31-mer probes
# (defined in 6c: present in the family arrays, absent from the other two
# families) directly in each library's reads via KMC.
#
# Outputs:
#   results/kmer_family_counts.csv   per family x sample: probe 31-mer occurrences,
#                                    normalized per 1e6 reads (probe_rpm)
#   results/kmer_family_fold.csv     CENP-A / H3K27ac fold + genomic-abundance ref
#
# Notes:
#   * probe sets are capped at 5000 31-mers/family (6c), so absolute RPM is
#     probe-set-size dependent; the CENP-A / control ratio is the meaningful
#     comparison, and the genomic-abundance line gives the no-enrichment null.
#   * no IgG/input exists in this dataset; H3K27ac is the negative control.
#
# Usage: Rscript 6e_kmer_analyze.R <work_dir>
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
PROBE <- file.path(WORK, "data", "probes")
RES   <- file.path(WORK, "results")
SRC   <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
dir.create(RES, showWarnings = FALSE)

SAMPLES <- c(XG_150="CENPA rep1", XG_151="CENPA rep2",
             XG_152="H3K27ac ctrl", XG_153="H3K27ac ctrl2")

# ---- library sizes (fragment BED line counts == read pairs) ----
lib_size <- c()
for (s in names(SAMPLES)) {
  frag <- file.path(SRC, "data", "fragments", paste0(s, "_fragments.bed"))
  lib_size[[s]] <- as.numeric(system(paste("wc -l <", frag), intern = TRUE))
}
cat("library sizes:", paste(names(lib_size), lib_size, sep="=", collapse=", "), "\n")

# ---- family bp (genomic abundance) ----
fam_bp <- c("349" = 113945325, "195" = 10781044, "389" = 2153924)  # chr-only TRF bin totals
genome_bp <- 3381500000   # chr1-28,X,Y total ~3.38 Gb
genomic_frac <- fam_bp / genome_bp

rows <- list()
for (fam in c("349", "195", "389")) {
  nprobes <- as.numeric(system(paste("wc -l <", file.path(PROBE, paste0("probes_", fam, ".31mers.txt"))), intern = TRUE))
  for (s in names(SAMPLES)) {
    f <- file.path(PROBE, paste0("probes_", fam, "_", s, "_counts.tsv"))
    if (!file.exists(f)) next
    d <- fread(f, header = FALSE, col.names = c("mer", "count"))
    probe_reads <- sum(d$count)
    probe_rpm   <- probe_reads / lib_size[[s]] * 1e6
    rows[[length(rows)+1]] <- data.table(
      family = fam, sample = s, sample_label = SAMPLES[[s]],
      n_probes = nprobes, n_probe_hits = nrow(d),
      probe_reads = probe_reads, probe_rpm = probe_rpm,
      lib_size = lib_size[[s]],
      genomic_frac = genomic_frac[[fam]],
      exp_reads_if_genomic = genomic_frac[[fam]] * lib_size[[s]])
  }
}
res <- rbindlist(rows)
fwrite(res, file.path(RES, "kmer_family_counts.csv"))

# ---- CENP-A vs control fold + genomic-abundance reference ----
wide <- dcast(res, family ~ sample, value.var = "probe_rpm")
wide[, cenpa_mean := (XG_150 + XG_151) / 2]
wide[, ctrl_mean  := (XG_152 + XG_153) / 2]
wide[, fold_vs_ctrl := cenpa_mean / ifelse(ctrl_mean == 0, NA, ctrl_mean)]
# expected probe RPM at genomic abundance, scaled by probe set coverage of the
# family is not directly estimable from RPM; report genomic_frac as the null
# "fraction of reads that should carry family sequence under no enrichment".
wide[, genomic_frac := genomic_frac[family]]
fwrite(wide, file.path(RES, "kmer_family_fold.csv"))

cat("=== kmer results (probe 31-mer occurrences per 1e6 reads) ===\n")
print(wide)
cat("DONE 6e_kmer_analyze.R\n")
