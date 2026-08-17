#!/usr/bin/env Rscript
# ============================================================================
# 195bp_l1_overlap.R — PLOTTING for the supplementary figure
#   "The 195 bp tandem repeats are L1 retrotransposon sequence"
#
# 3 panels (a-b top, c full width bottom):
#   (a) RepeatMasker class composition of the bin4 (195 bp) vs bin6 (349 bp)
#       arrays — an EXCLUSIVE partition (L1 > Unknown > Other > Unannotated),
#       so each bar sums to 100% of array bp. Shows 195 bp = 99.4% LINE/L1
#       (fam-189/18) vs 349 bp = 99.3% Unknown (species-specific satellite).
#   (b) The 195 bp unit is the L1's 5'-terminal tandem repeat: rnd-1_family-189
#       and rnd-1_family-18 consensus with the 195 bp repeat blocks at the
#       5' end (verified by TRF, same params as the genome call).
#   (c) Example locus chr21:83.26-83.34 Mb — the longest 195 bp array is a
#       tandem array of L1 fragments (repeating unit: fam-189 + fam-18 +
#       LTR/other + Unknown junctions), with the TRF 195 bp array below.
#
# Reads cached tables from results/ (built by prepare_panel_data.py).
# Outputs 195bp_l1_overlap.{pdf,png,svg}  (svg via svglite for Illustrator).
# FONT_SIZE is the single manual tuning knob.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

FONT_SIZE <- 12
OUT_PREFIX <- "195bp_l1_overlap"
OUT_DIR <- "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/195bp-repeatmasker-overlap"
RES_DIR  <- file.path(OUT_DIR, "results")

# ── Read cached data ────────────────────────────────────────────────────────
comp <- fread(file.path(RES_DIR, "panel_a_composition.tsv"))
cons <- fread(file.path(RES_DIR, "panel_b_consensus.tsv"))
locus <- fread(file.path(RES_DIR, "panel_c_locus.tsv"))
arrs  <- fread(file.path(RES_DIR, "panel_c_arrays.tsv"))

# ── Shared category palette ─────────────────────────────────────────────────
CAT_LEVELS <- c("l1_fam189", "l1_fam18", "l1_other", "unknown", "other", "unannotated")
CAT_LABELS <- c("L1 rnd-1_family-189", "L1 rnd-1_family-18", "L1 other",
                "Unknown (no RM match)", "Other (LTR/SINE/DNA/simple)", "Unannotated")
CAT_COLORS <- c(l1_fam189 = "#2166AC", l1_fam18 = "#74ADD1", l1_other = "#B3CDE3",
                unknown = "#8C8C8C", other = "#D8C9A3", unannotated = "#F2F2F2")

# ════════════════════════════════════════════════════════════════════════════
# PANEL (a) — class composition
# ════════════════════════════════════════════════════════════════════════════
comp[, category := factor(category, levels = CAT_LEVELS)]
comp[, array_group := factor(array_group, levels = c("bin4_195bp", "bin6_349bp"))]
grp_lab <- c(bin4_195bp = "195 bp arrays\n(12,631 arrays · 10.6 Mb)",
             bin6_349bp = "349 bp arrays\n(799 arrays · 113.3 Mb)")
# L1 % per group for the annotation
l1pct <- comp[category %in% c("l1_fam189", "l1_fam18", "l1_other"),
              .(l1 = sum(pct)), by = array_group]
ann <- data.table(
  array_group = c("bin4_195bp", "bin6_349bp"),
  lab = c("99.4% LINE/L1", "99.3% no RepeatMasker match\n(species-specific satellite)")
)

pA <- ggplot(comp, aes(x = array_group, y = pct, fill = category)) +
  geom_bar(stat = "identity", width = 0.62, colour = NA) +
  geom_text(data = ann, aes(x = array_group, y = 100, label = lab),
            inherit.aes = FALSE, hjust = 1.02, size = 0.30 * FONT_SIZE,
            fontface = "plain", lineheight = 0.9) +
  scale_fill_manual(values = CAT_COLORS, labels = CAT_LABELS, drop = FALSE) +
  scale_x_discrete(labels = grp_lab) +
  scale_y_continuous(limits = c(0, 103), expand = c(0, 0),
                     labels = function(x) paste0(x, "%")) +
  coord_flip() +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_minimal(base_size = FONT_SIZE) +
  theme(
    axis.text.y = element_text(size = FONT_SIZE - 0.5, color = "black"),
    axis.text.x = element_text(size = FONT_SIZE - 1.5, color = "black"),
    axis.title  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    legend.position = "bottom",
    legend.text  = element_text(size = FONT_SIZE - 1.5),
    legend.title = element_blank(),
    legend.key.size = unit(0.32, "cm"),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(6, 6, 2, 6)
  )

# ════════════════════════════════════════════════════════════════════════════
# PANEL (b) — L1 consensus 5' tandem repeat
# ════════════════════════════════════════════════════════════════════════════
# Parse the unit blocks into segments
blk <- list()
for (i in seq_len(nrow(cons))) {
  fam <- cons$family[i]
  for (b in strsplit(cons$blocks_csv[i], "\\|")[[1]]) {
    ss <- as.integer(strsplit(b, "-")[[1]])
    blk[[length(blk) + 1]] <- data.table(family = fam, start = ss[1], end = ss[2])
  }
}
blk <- rbindlist(blk)
blk[, y := ifelse(family == "rnd-1_family-189", 2L, 1L)]
blk[, shade := rep(c("#2166AC", "#74ADD1"), length.out = .N), by = family]
# count FULL (195 bp) copies per family for labels
blk[, is_partial := (end - start + 1) < 195]
nfull <- blk[is_partial == FALSE, .(n = .N), by = family]
cons <- merge(cons, nfull, by = "family")
cons[, y := ifelse(family == "rnd-1_family-189", 2L, 1L)]
cons[, lab := paste0(family, "  ·  ", consensus_len, " bp  ·  ",
                     n, " × 195 bp + partial")]

pB <- ggplot() +
  # full consensus body (light grey)
  geom_segment(data = cons,
               aes(x = 0, xend = consensus_len, y = y, yend = y),
               linewidth = 8, colour = "#DCDCDC", lineend = "butt") +
  # 5' repeat blocks
  geom_rect(data = blk,
            aes(xmin = start, xmax = end, ymin = y - 0.42, ymax = y + 0.42,
                fill = shade), colour = "white", linewidth = 0.3) +
  scale_fill_identity() +
  # 5' / 3' labels
  annotate("text", x = -700, y = 2, label = "5′", hjust = 1,
           size = 0.30 * FONT_SIZE, fontface = "italic", colour = "grey25") +
  annotate("text", x = -700, y = 1, label = "5′", hjust = 1,
           size = 0.30 * FONT_SIZE, fontface = "italic", colour = "grey25") +
  annotate("text", x = max(cons$consensus_len) + 700, y = 2, label = "3′",
           hjust = 0, size = 0.30 * FONT_SIZE, fontface = "italic", colour = "grey25") +
  annotate("text", x = max(cons$consensus_len) + 700, y = 1, label = "3′",
           hjust = 0, size = 0.30 * FONT_SIZE, fontface = "italic", colour = "grey25") +
  # family labels (right side)
  geom_text(data = cons, aes(x = consensus_len + 1250, y = y, label = lab),
            hjust = 0, size = 0.30 * FONT_SIZE, colour = "black") +
  # 5' repeat bracket + caption
  annotate("segment", x = 20, xend = blk[family == "rnd-1_family-189", max(end)],
           y = 2.78, yend = 2.78, linewidth = 0.4, colour = "grey30") +
  annotate("text", x = (20 + blk[family == "rnd-1_family-189", max(end)]) / 2,
           y = 2.92, label = "5′ tandem repeat  (195 bp unit)",
           size = 0.28 * FONT_SIZE, colour = "grey30") +
  scale_x_continuous(limits = c(-1400, max(cons$consensus_len) + 6200),
                     expand = expansion(mult = c(0, 0.01)),
                     breaks = seq(0, 8000, 2000)) +
  scale_y_continuous(limits = c(0.35, 3.2), expand = c(0, 0)) +
  labs(x = "L1 consensus position (bp)") +
  theme_minimal(base_size = FONT_SIZE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = FONT_SIZE - 2, color = "black"),
    axis.title.x = element_text(size = FONT_SIZE - 1, margin = margin(t = 2), color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey85", linewidth = 0.3),
    plot.margin = margin(6, 6, 2, 6)
  )

# ════════════════════════════════════════════════════════════════════════════
# PANEL (c) — example locus track
# ════════════════════════════════════════════════════════════════════════════
locus[, category := factor(category, levels = c("L1 fam-189", "L1 fam-18", "L1 other",
                                                "Unknown", "Other"))]
locus[, y := 2]
arrs[, y := 1]
LOC_COLORS <- c("L1 fam-189" = "#2166AC", "L1 fam-18" = "#74ADD1", "L1 other" = "#B3CDE3",
                "Unknown" = "#8C8C8C", "Other" = "#D8C9A3")
xmin <- 83255000; xmax <- 83344000
# give room on the left for the track labels
xlims <- c(xmin - 60000, xmax + 5000)

pC <- ggplot() +
  geom_rect(data = locus,
            aes(xmin = start, xmax = end, ymin = y - 0.30, ymax = y + 0.30, fill = category),
            colour = NA) +
  geom_rect(data = arrs,
            aes(xmin = start, xmax = end, ymin = y - 0.42, ymax = y + 0.42),
            fill = "#111111", colour = NA) +
  scale_fill_manual(values = LOC_COLORS, drop = FALSE) +
  # track labels (left margin)
  annotate("text", x = xmin - 15000, y = 2, label = "RepeatMasker", hjust = 1,
           size = 0.30 * FONT_SIZE, colour = "black") +
  annotate("text", x = xmin - 15000, y = 1, label = "195 bp\nTRF arrays", hjust = 1,
           size = 0.30 * FONT_SIZE, colour = "black", lineheight = 0.9) +
  # scale bar
  annotate("segment", x = xmax - 50000, xend = xmax - 50000 + 10000, y = 0.25, yend = 0.25,
           linewidth = 0.5, colour = "black") +
  annotate("text", x = xmax - 50000 + 5000, y = 0.12, label = "10 kb",
           size = 0.26 * FONT_SIZE, colour = "black") +
  scale_x_continuous(limits = xlims, expand = c(0, 0),
                     labels = function(x) paste0(sprintf("%.3f", x / 1e6), " Mb"),
                     breaks = seq(83.26e6, 83.44e6, by = 0.04e6)) +
  scale_y_continuous(limits = c(0.05, 2.75), expand = c(0, 0)) +
  labs(x = "chr21 position") +
  theme_minimal(base_size = FONT_SIZE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = FONT_SIZE - 2, color = "black", angle = 0, hjust = 0.5),
    axis.title.x = element_text(size = FONT_SIZE - 1, margin = margin(t = 2), color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey85", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(6, 6, 2, 6)
  )

# ── Compose + save ──────────────────────────────────────────────────────────
fig <- (pA | pB) / pC + plot_layout(heights = c(1, 1))

ggsave(file.path(OUT_DIR, paste0(OUT_PREFIX, ".pdf")), fig,
       width = 11, height = 6.2, device = cairo_pdf)
ggsave(file.path(OUT_DIR, paste0(OUT_PREFIX, ".png")), fig,
       width = 11, height = 6.2, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, paste0(OUT_PREFIX, ".svg")), fig,
       width = 11, height = 6.2, dpi = 300, bg = "white", device = svglite::svglite)

cat("wrote", file.path(OUT_DIR, OUT_PREFIX), "{pdf,png,svg}\n")
