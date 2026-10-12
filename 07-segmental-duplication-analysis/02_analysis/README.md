# 07-segmental-duplication-analysis

Detect segmental duplications (SDs) with BISER and produce the Figure 5 analysis
(cumulative-length "pervasive" region, gene overlap, GO enrichment, pyCirclize
circos, and paralog expression).

## Order of operation

| Step | File | What it does |
|---|---|---|
| 1 | `01_biser_run.sh` | run BISER on one genome → `segdup_output_<species>.bedpe` |
| 2 | `02_biser_submit_batch.sh` | submit step 1 for every genome in `BISER_GENOME_LIST` (optional) |
| 2b | `02_segdup_per_chromosome.ipynb` | per-chromosome segdup summary (intra/inter counts, scores, merged coverage, % of chromosome) |
| 3 | `03_understanding_segdup_output.ipynb` | inspect raw BEDPE, drop reciprocal duplicates → `segdup_output_duplicateLinkRemoved.bedpe` |
| 4 | `04_segdup_gene_overlap.ipynb` | gene × SD overlap (agat GFF) → `hifiasm_gene_segDup_overlapInfo_092525.tsv` |
| 5 | `05_cumulative_feature_overlap.ipynb` | cumulative-length plot, 80% "pervasive" region, feature overlap (input to **Fig 5A**) |
| 6 | `06_go_enrichment.ipynb` | GO enrichment of pervasive-region genes (input to **Fig 5B**) |
| 7 | `07_separate_gff.ipynb` | split the agat GFF into `nc` / `gr` / `ga` gene tracks |
| 8 | `08_circos_pyCirclize.ipynb` | pyCirclize circos of pervasive region + genes of interest (earlier circos; final panel: `11b`) |
| 9 | `09_rna_seq_level.sh` | STAR index / RNA-seq level inputs (CCT7 paralogs) |
| 10 | `10_expression_heatmap_R.ipynb` | paralog expression heatmap across 29 samples / 5 organs (earlier heatmap version; not in the final figure) |

### Figure 5 panel scripts (printed size, 8 pt; outputs in `../figure5/`)

| Panel | Script | Reads |
|---|---|---|
| A | `05b_plot_figure5A_cumulative.py` | `merged_segdup_with_feature_counts.tsv` (step 5 cache) + `$SEGDUP_GENE_OVERLAP_TSV` |
| B | `06b_plot_figure5B_go_dotplot.py` | `go_enrichment_results_segdup_hotspot_seqRemoved_80percLength.csv` (step 6) |
| C | `11b_plot_figure5C_circos.py` | BISER BEDPE, `segdup_density_rolling_5000000_100000_bpLevel_hotspot80perc.tsv`, `segdup_genes_hotspots_80percLength.tsv`, `merged_paralog_singlecopy_families.tsv` |
| D | `12b_plot_figure5D_expression_violin.py` | `merged_paralog_singlecopy_counts.tsv`, `merged_paralog_singlecopy_families.tsv` |

"Largest SegDups" (file names say `hotspot80perc` / `80percLength`) = the 1,425 longest merged
regions whose cumulative length is <= 80% of total SegDup length (28,961 regions); they overlap
1,271 genes (469 non-LOC symbols used for GO enrichment). No evolutionary "hotspot" claim is implied.

## Inputs

- **Segdup BEDPE**: `.../biser/hifiasm-041425/segdup_output.bedpe` (step 3) and
  `.../segdup_output_duplicateLinkRemoved.bedpe` (steps 4–5); see `BISER_BEDPE`
  in `config.sh`.
- **Gene annotation** (the new merged annotation, AGAT-normalized):
  `code/command-line-script/annotation-merging/output/hifiasm-041425-denovoEnhanced_parentReassigned_20260820.gff3`
  — used by steps 4 and 7. Gene IDs are `ID=gene-…`.
- **Repeat annotation** (RepeatMasker `.out.gff`) for the feature-overlap in step 5.

Steps 4–5 must be re-run after a new annotation (step 6 consumes step 4's output;
step 7 must also be re-run to refresh the `nc`/`gr`/`ga` tracks for step 8).

## Configuration

BISER run settings (`BISER_BIN`, `BISER_OUT_DIR`, `BISER_THREADS`, `BISER_GC_HEAP`,
`BISER_MAX_EDIT_ERROR`, `BISER_MAX_ERROR`, `BISER_GENOME_DIR`, `BISER_GENOME_LIST`)
are in `config.sh`. The `genome-annotation` conda env (`ENV_GENOME_ANNOTATION`)
supplies `samtools` for indexing.

## Notes

- BISER (v1.4) is a user-level Python install: entry point `~/.local/bin/biser`
  (`BISER_BIN`), package in `~/.local/lib/python3.9/site-packages/biser`.
- Step 7 reads three gene lists in
  `code/jupyter-notebook-script/gff-build/final-list/` (`gene_LOC_nameChange_unique.tsv`
  → `nc`, `gene_LOC_replace_unique.tsv` → `gr`, `gene_gBraker_toAdd_unique.tsv` → `ga`);
  these were generated from the *old* annotation and may need regenerating after a
  re-merge.
- The notebooks were cleaned of exploratory cells and have their old outputs
  cleared; re-run them in the numbered order after updating inputs.
- `02_segdup_per_chromosome.ipynb` reads `segdup_output_duplicateLinkRemoved.bedpe`
  (the output of step 3). The archived source it was pulled from was missing two
  cells — the `pivot_lengths` pivot and the per-chromosome `result_df` totals —
  which were reconstructed here.
- `#SBATCH` scripts target TSCC/UCSD; adjust for your cluster.
