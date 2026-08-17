#!/usr/bin/env Rscript
# ============================================================================
# 8_plots.R -- figures for the CENP-A-domain-centered repeat-composition analysis
#
# Main:
#   Panel A: 9-bin stacked repeat composition by distance band, per-chromosome
#            (spaghetti lines, one line per chromosome, y = band fraction)
#   Panel B: 9-bin stacked composition, genomewide pooled (all chromosomes
#            combined into one CENP-A peak) -- stacked bars per band
#   Panel C: family concentration (frac of family bp in each band) observed vs
#            matched-null, per family (349/195/389)
#
# Supplementary:
#   Supp A: domain definition comparison (weighted vs k1), esp. chr4
#   Supp B: matched-null covariate QC
#   Supp C: k-mer family content
#   Supp D: unique vs weighted CENP-A signal by band
#
# Usage: Rscript 8_plots.R <work_dir>
# ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
WORK <- if (length(args) >= 1) args[1] else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
DATA <- file.path(WORK, "data")
RES  <- file.path(WORK, "results")
MAIN <- file.path(WORK, "plots", "main")
SUPP <- file.path(WORK, "plots", "supp")
dir.create(MAIN, showWarnings = FALSE); dir.create(SUPP, showWarnings = FALSE)

BAND_ORDER <- c("core", "core_to_250kb", "250kb_to_1Mb", "remainder")
BAND_LAB   <- c(core = "CENP-A core", core_to_250kb = "core→250 kb",
                `250kb_to_1Mb` = "250 kb→1 Mb", remainder = "remainder")
BIN9_COLS <- c("#8c510a","#d8b365","#f6e8c3","#c7eae5","#5ab4ac",
               "#2166ac","#762a83","#e7298a","#1b7837")
names(BIN9_COLS) <- c("1-10 bp","11-50 bp","51-192 bp","193-195 bp","196-347 bp",
                      "348-349 bp","350-385 bp","386-390 bp","391+ bp")
BIN9_COLS["other (no TRF)"] <- "#d9d9d9"

# ---- load per-chromosome 9-bin composition ----
comp <- fread(file.path(RES, "repeat_composition_by_band_9bin.csv"))
comp[, band := factor(band, levels = BAND_ORDER)]

# ---- Panel A: per-chromosome spaghetti, 9 bins (facets) ----
pA <- ggplot(comp, aes(x = band, y = frac_band * 100, group = chrom, color = chrom)) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(~ bin_label, ncol = 5) +
  scale_x_discrete(labels = function(x) sub("core→", "", x)) +
  labs(title = "9-bin TRF composition by distance from CENP-A domain (per chromosome)",
       x = "distance band", y = "% of band bp covered",
       color = "chromosome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")
ggsave(file.path(MAIN, "panelA_9bin_per_chromosome.pdf"), pA, width = 11, height = 6)

# ---- Panel B: genomewide pooled stacked bars ----
gw <- fread(file.path(RES, "repeat_composition_genomewide_9bin.csv"))
gw[, band := factor(band, levels = BAND_ORDER)]
# keep order by increasing bin for a sensible stack
gw[, bin_order := ifelse(bin_id == 0, 10, bin_id)]
pB <- ggplot(gw, aes(x = band, y = frac_band * 100, fill = bin_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = BIN9_COLS) +
  scale_x_discrete(labels = BAND_LAB[levels(gw$band)]) +
  labs(title = "Genomewide 9-bin repeat composition by distance from CENP-A peak\n(all chromosomes pooled into one CENP-A domain)",
       x = "distance band", y = "% of band bp", fill = "TRF period") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(MAIN, "panelB_9bin_genomewide_pooled.pdf"), pB, width = 8, height = 6)

# ---- Panel C: family concentration observed vs matched null ----
enr <- fread(file.path(RES, "matched_null_band_enrichment.csv"))
enr[, band := factor(band, levels = BAND_ORDER)]
enr_long <- melt(enr[, .(family, band, obs_mean, null_mean)],
                 id.vars = c("family", "band"),
                 measure.vars = c("obs_mean", "null_mean"),
                 variable.name = "type", value.name = "frac")
enr_long[type == "obs_mean", type := "observed"]
enr_long[type == "null_mean", type := "matched-null"]
pC <- ggplot(enr_long, aes(x = band, y = frac * 100, fill = type, group = type)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  facet_wrap(~ family) +
  scale_x_discrete(labels = BAND_LAB[levels(enr_long$band)]) +
  labs(title = "Family bp fraction in each distance band: observed vs matched-shuffle null",
       x = "distance band", y = "% of family bp", fill = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")
ggsave(file.path(MAIN, "panelC_family_vs_null.pdf"), pC, width = 10, height = 5)

# ---- Supp A: domain definition comparison ----
for (sig in c("weighted", "k1")) {
  f <- file.path(DATA, "domains", paste0("cenpa_domains_", sig, ".csv"))
  if (!file.exists(f)) next
  d <- fread(f)
  d[, sig := sig]
  assign(paste0("dom_", sig), d)
}
if (exists("dom_weighted") && exists("dom_k1")) {
  dd <- rbind(dom_weighted, dom_k1)
  pS1 <- ggplot(dd, aes(x = chrom, y = core_size / 1e6, fill = sig)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    labs(title = "CENP-A core size per chromosome (weighted k=100 vs k=1)",
         x = "chromosome", y = "core size (Mb)", fill = "domain basis") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(file.path(SUPP, "suppA_domain_weighted_vs_k1.pdf"), pS1, width = 12, height = 4)
}

# ---- Supp C: k-mer family content ----
kf <- fread(file.path(RES, "kmer_family_fold.csv"))
pS3 <- ggplot(kf, aes(x = family, y = cenpa_mean, fill = "CENP-A")) +
  geom_col() +
  geom_point(aes(y = ctrl_mean), shape = 21, size = 3, fill = "grey70") +
  labs(title = "Family-specific 31-mer content in CENP-A fragments\n(occurrences per 1e6 reads; dots = H3K27ac control)",
       x = "repeat family", y = "probe 31-mer occurrences / 1e6 reads") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(SUPP, "suppC_kmer_family_content.pdf"), pS3, width = 5, height = 4)

# ---- Supp D: unique vs weighted signal by band ----
sig <- fread(file.path(RES, "cenpa_signal_by_band_window.csv"))
sig[, band := factor(band, levels = BAND_ORDER)]
sig_long <- melt(sig[, .(chrom, band, weighted_mean, k1_mean)],
                 id.vars = c("chrom", "band"),
                 measure.vars = c("weighted_mean", "k1_mean"),
                 variable.name = "type", value.name = "signal")
pS4 <- ggplot(sig_long, aes(x = band, y = signal, fill = type)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_x_discrete(labels = BAND_LAB[levels(sig_long$band)]) +
  labs(title = "CENP-A signal by distance band: unique (k=1) vs 1/NH-weighted k=100",
       x = "distance band", y = "mean CENP-A signal (CPM)", fill = "mapping") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(SUPP, "suppD_unique_vs_weighted.pdf"), pS4, width = 8, height = 5)

cat("DONE 8_plots.R\n")
