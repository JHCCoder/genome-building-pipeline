#!/usr/bin/env Rscript
# Figure 2D/E, step 1: cluster the enriched GO BP terms and cache everything both panels need.
#
# Same steps as 04_go_pcoa_barplot_R.ipynb (cells 3-9, 13): significant BP terms
# (adj. P <= 0.05) from the degu-background enrichment -> rrvgo "Rel" similarity ->
# reduceSimMatrix(threshold = 0.5, scores = -log10 adj. P) -> 2-D PCoA (cmdscale on
# 1 - similarity). The top 9 clusters are those whose most significant member term has the
# highest -log10(adj. P) (chosen 2026-10-11; the notebook's "top 5 by overlap fraction" tied
# six clusters at 1.0). D colours these clusters; E colours its bars by the same clusters.
# Each top cluster's colour is stored in the cache so both panels use the same one.
#
# Run with the r-clusterProfiler env:
#   Rscript 04_prepare_figure2DE_go_clusters.R

suppressPackageStartupMessages({
  library(rrvgo)
  library(org.Hs.eg.db)
  library(dplyr)
})

.cfg <- function(name) {
  v <- Sys.getenv(name)
  if (nzchar(v)) return(v)
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  root <- normalizePath(if (length(arg)) dirname(sub("^--file=", "", arg[1])) else getwd())
  while (!file.exists(file.path(root, "config.sh")) && root != dirname(root)) root <- dirname(root)
  cmd <- sprintf("source \"%s\" 2>/dev/null; printf \"%%s\" \"${%s}\"",
                 file.path(root, "config.sh"), name)
  v <- system2("bash", c("-c", shQuote(cmd)), stdout = TRUE)
  if (length(v) == 0 || !nzchar(v)) stop(sprintf("ERROR: %s is not set in the environment or in config.sh", name))
  v
}

GO_TYPE   <- "human_GO_bp_2025.gmt"
THRESHOLD <- 0.05
N_TOP     <- 9
# Paul Tol "muted": nine colour-blind-safe categorical colours (Okabe-Ito has only eight)
PALETTE   <- c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", "#999933", "#AA4499")

go_df <- read.csv(.cfg("FIG2_GO_TSV"), sep = "\t", stringsAsFactors = FALSE)[, -1]
go_int <- go_df[go_df$Gene_set == GO_TYPE & go_df$Adjusted.P.value <= THRESHOLD, ]
cat("significant BP terms:", nrow(go_int), "\n")

sim <- calculateSimMatrix(go_int$GO_ID, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
scores <- setNames(-log10(go_int$Adjusted.P.value), go_int$GO_ID)
reduced <- reduceSimMatrix(sim, scores = scores, threshold = 0.5, orgdb = "org.Hs.eg.db")

pcoa <- cmdscale(as.matrix(as.dist(1 - sim)), eig = TRUE, k = 2)$points
coords <- data.frame(go = rownames(pcoa), PCoA1 = pcoa[, 1], PCoA2 = pcoa[, 2])

out <- reduced %>%
  select(go, term, parent, parentTerm, score, cluster) %>%
  inner_join(coords, by = "go") %>%
  left_join(go_int %>% select(GO_ID, Overlap, Adjusted.P.value), by = c("go" = "GO_ID")) %>%
  mutate(overlap_fraction = sapply(strsplit(Overlap, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
         term = gsub("\\s*\\(GO:\\d+\\)", "", term),
         parentTerm = gsub("\\s*\\(GO:\\d+\\)", "", parentTerm))

top <- out %>% group_by(parent, parentTerm) %>%
  summarise(best_score = max(score), n_terms = n(), .groups = "drop") %>%
  arrange(desc(best_score)) %>% head(N_TOP) %>%
  mutate(top_rank = row_number(), color = PALETTE[top_rank])
out <- out %>% left_join(top %>% select(parent, top_rank, color), by = "parent")

dir.create(dirname(.cfg("FIG2DE_CLUSTERS_TSV")), showWarnings = FALSE, recursive = TRUE)
write.table(out, .cfg("FIG2DE_CLUSTERS_TSV"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("terms in PCoA:", nrow(out), " clusters:", length(unique(out$parent)), "\n")
cat("top", N_TOP, "clusters (by best -log10 adj. P):\n")
print(as.data.frame(top))
cat("terms in the top clusters:", sum(!is.na(out$top_rank)), "\n")
