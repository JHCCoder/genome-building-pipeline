#!/usr/bin/env Rscript
# ===========================================================================
# Heatmap of paralog expression (TPM) across 29 RNA-seq samples.
# Reimplementation using ComplexHeatmap (R/Bioconductor).
#
# - Rows: 2,252 paralogs, clustered within family blocks
# - Columns: 29 samples, grouped by tissue
# - Right annotation: top-40 family labels with connecting lines (anno_mark)
# - Top annotation: tissue color bar + tissue names
# ===========================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
ht_opt$message <- FALSE

# ===========================================================================
# Load data
# ===========================================================================
counts_raw <- read.table('paralog_families_counts.tsv', sep = '\t',
                         header = TRUE, comment.char = '#',
                         row.names = 1, check.names = FALSE)

# Extract count columns (SRR*)
sample_cols_orig <- grep('^SRR', colnames(counts_raw), value = TRUE)
counts <- counts_raw[, sample_cols_orig, drop = FALSE]
colnames(counts) <- gsub('_Aligned.sortedByCoord.out.bam', '', colnames(counts))
sample_list <- colnames(counts)

# Load family metadata
families_df <- read.table('paralog_families.tsv', sep = '\t', header = TRUE,
                          stringsAsFactors = FALSE)
rownames(families_df) <- families_df$gene_id

# Merge and filter to paralogs only
meta_cols <- c('family', 'gene_name', 'gene_type', 'paralog_type', 'length')
merged <- merge(counts, families_df[, meta_cols, drop = FALSE], by = 'row.names')
rownames(merged) <- merged$Row.names
merged$Row.names <- NULL

paralogs <- merged[merged$gene_type != 'parent', ]
cat(sprintf('Paralogs: %d\n', nrow(paralogs)))

# ===========================================================================
# TPM normalization
# ===========================================================================
gene_lengths <- paralogs$length
rpk <- paralogs[, sample_list] / (gene_lengths / 1000)
rpk <- as.matrix(rpk)
tpm <- sweep(rpk, 2, colSums(rpk), '/') * 1e6

# ===========================================================================
# Sample / tissue metadata
# ===========================================================================
sample_tissue <- c(
  'SRR17216293' = 'skin',    'SRR17216294' = 'brain',   'SRR17216295' = 'heart',
  'SRR17216296' = 'kidney',  'SRR17216297' = 'liver',   'SRR17216298' = 'lung',
  'SRR17216299' = 'skin',    'SRR17216300' = 'brain',   'SRR17216301' = 'heart',
  'SRR17216302' = 'kidney',  'SRR17216303' = 'liver',   'SRR17216304' = 'lung',
  'SRR17216305' = 'skin',    'SRR17216306' = 'brain',   'SRR17216307' = 'heart',
  'SRR17216308' = 'kidney',  'SRR17216309' = 'liver',   'SRR17216310' = 'lung',
  'SRR17216311' = 'brain',   'SRR17216312' = 'heart',   'SRR17216313' = 'kidney',
  'SRR17216314' = 'kidney',  'SRR17216315' = 'lung',    'SRR17216316' = 'skin',
  'SRR17216317' = 'heart',   'SRR17216318' = 'kidney',  'SRR17216319' = 'liver',
  'SRR17216320' = 'brain',   'SRR17216321' = 'lung'
)
tissues <- c('skin', 'brain', 'heart', 'kidney', 'liver', 'lung')
tissue_colors <- c(
  'skin'   = '#E8A87C',
  'brain'  = '#C38D9E',
  'heart'  = '#E27D60',
  'kidney' = '#85CDCA',
  'liver'  = '#D4A574',
  'lung'   = '#41B3A3'
)

# Sort samples by tissue
tissue_order <- match(sample_tissue[sample_list], tissues)
sample_order_idx <- order(tissue_order)
sample_order <- sample_list[sample_order_idx]
tpm_sorted <- tpm[, sample_order]

# ===========================================================================
# Row ordering: cluster families, keep family members together
# ===========================================================================
# Compute per-family mean expression vector across samples
family_mean <- aggregate(tpm, by = list(family = paralogs$family), FUN = mean)
rownames(family_mean) <- family_mean$family
family_mean$family <- NULL
family_mean <- as.matrix(family_mean)

# Cluster families (ward linkage, matches Python scipy linkage(method='ward'))
family_dist <- dist(family_mean)
family_hclust <- hclust(family_dist, method = 'ward.D2')
family_order <- family_hclust$labels[family_hclust$order]

# Within each family, sort paralogs by total expression (descending)
row_order <- c()
family_boundaries <- list()

for (fam in family_order) {
  fam_mask <- paralogs$family == fam
  fam_genes <- rownames(paralogs)[fam_mask]
  fam_expr <- rowMeans(tpm_sorted[fam_genes, , drop = FALSE])
  fam_order_genes <- names(sort(fam_expr, decreasing = TRUE))

  start <- length(row_order) + 1
  row_order <- c(row_order, fam_order_genes)
  end <- length(row_order)
  family_boundaries[[fam]] <- c(start = start, end = end)
}

tpm_ordered <- tpm_sorted[row_order, ]

# Log-transform for visualization
log_tpm <- log2(as.matrix(tpm_ordered) + 0.1)

cat(sprintf('Heatmap dimensions: %d rows x %d columns\n',
            nrow(log_tpm), ncol(log_tpm)))
cat(sprintf('Families: %d\n', length(family_order)))

# ===========================================================================
# Top-40 families for right-side annotation
# ===========================================================================
family_sizes <- table(paralogs$family)
top40_fams <- names(sort(family_sizes, decreasing = TRUE)[1:min(40, length(family_sizes))])

# Build anno_mark data: midpoints + labels for each top-40 family
mark_at <- c()
mark_labels <- c()

for (fam in family_order) {
  if (fam %in% top40_fams) {
    bounds <- family_boundaries[[fam]]
    mid <- (bounds['start'] + bounds['end']) / 2
    n_members <- bounds['end'] - bounds['start'] + 1
    mark_at <- c(mark_at, mid)
    mark_labels <- c(mark_labels, sprintf('%s  (%d)', fam, n_members))
  }
}

cat(sprintf('Annotated families: %d\n', length(mark_at)))

# ===========================================================================
# Column annotation (tissue bar + tissue names)
# ===========================================================================
col_tissue <- sample_tissue[sample_order]
col_tissue_factor <- factor(col_tissue, levels = tissues)

# Compute replicate count label per tissue
tissue_counts <- table(col_tissue_factor)
rep_labels <- sprintf('%d rep.', tissue_counts[tissues])

# Top annotation: colored bar for each tissue
col_anno <- HeatmapAnnotation(
  Tissue = col_tissue_factor,
  col = list(Tissue = tissue_colors),
  show_legend = FALSE,
  annotation_name_side = 'left',
  annotation_name_gp = gpar(fontsize = 0),  # hide annotation name
  simple_anno_size = unit(0.35, 'cm'),
  border = FALSE,
  gap = unit(0, 'mm')
)

# ===========================================================================
# Color mapping (YlOrRd-like)
# ===========================================================================
vmax <- quantile(log_tpm, 0.99)
col_fun <- colorRamp2(
  seq(-3, vmax, length.out = 9),
  c('#FFFFCC', '#FFEDA0', '#FED976', '#FEB24C', '#FD8D3C',
    '#FC4E2A', '#E31A1C', '#BD0026', '#800026')
)

# ===========================================================================
# Right annotation: family labels with connecting lines
# ===========================================================================
right_anno <- rowAnnotation(
  Family = anno_mark(
    at = mark_at,
    labels = mark_labels,
    which = 'row',
    side = 'right',
    labels_gp = gpar(fontsize = 9, fontface = 'bold', col = '#222222'),
    link_width = unit(4, 'cm'),
    link_gp = gpar(col = '#AAAAAA', lwd = 0.5),
    extend = unit(0, 'mm')
  )
)

# ===========================================================================
# Main heatmap
# ===========================================================================
ht <- Heatmap(
  log_tpm,
  name = 'log2TPM',   # simple name — no special chars (legend title set separately)

  # Color
  col = col_fun,

  # No built-in clustering (we ordered manually)
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,

  # Column splits by tissue (adds white separators)
  column_split = col_tissue_factor,
  column_gap = unit(1.5, 'mm'),
  column_title = tissues,
  column_title_side = 'top',
  column_title_gp = gpar(fontsize = 14, fontface = 'bold', col = 'black'),

  # Annotations
  top_annotation = col_anno,
  right_annotation = right_anno,

  # Visual polish
  border = '#444444',
  rect_gp = gpar(col = NA),
  row_title = NULL,
  column_names_rot = 0,

  # Do not draw dendrograms
  show_row_dend = FALSE,
  show_column_dend = FALSE,

  # Heatmap legend
  heatmap_legend_param = list(
    title = expression(log[2] ~ (TPM + 0.1)),
    title_position = 'lefttop-rot',
    title_gp = gpar(fontsize = 11),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(4, 'cm'),
    direction = 'vertical'
  )
)

# ===========================================================================
# Draw and save
# ===========================================================================

# Build combined title with subtitle
title_str <- sprintf(
  'Paralog expression across tissues\n%s paralogs | %s families | %d RNA-seq samples | Top 40 families labeled',
  format(nrow(paralogs), big.mark = ','),
  format(length(family_order), big.mark = ','),
  ncol(log_tpm))

# Helper to draw the full figure
draw_full <- function() {
  draw(ht,
       merge_legend = TRUE,
       heatmap_legend_side = 'left',
       padding = unit(c(2, 10, 2, 2), 'cm'),
       column_title = title_str,
       column_title_gp = gpar(fontsize = 18, fontface = 'bold'))
}

# PDF
pdf('paralog_expression_heatmap_all.pdf', width = 16, height = 12)
draw_full()
dev.off()
cat('Saved: paralog_expression_heatmap_all.pdf\n')

# PNG
png('paralog_expression_heatmap_all.png', width = 16, height = 12,
    units = 'in', res = 600)
draw_full()
dev.off()
cat('Saved: paralog_expression_heatmap_all.png\n')

cat('Done.\n')
