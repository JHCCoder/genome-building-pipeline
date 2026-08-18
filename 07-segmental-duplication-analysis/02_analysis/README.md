# 06-segmental-duplication-analysis

Detect segmental duplications (SDs) with BISER and produce the Figure 7 analysis
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
| 5 | `05_cumulative_feature_overlap.ipynb` | cumulative-length plot, 80% "pervasive" region, feature overlap (**Fig 7A**) |
| 6 | `06_go_enrichment.ipynb` | GO enrichment of pervasive-region genes (**Fig 7B**) |
| 7 | `07_separate_gff.ipynb` | split the agat GFF into `nc` / `gr` / `ga` gene tracks |
| 8 | `08_circos_pyCirclize.ipynb` | pyCirclize circos of pervasive region + genes of interest (**Fig 7C**) |
| 9 | `09_rna_seq_level.sh` | STAR index / RNA-seq level inputs (CCT7 paralogs) |
| 10 | `10_expression_heatmap_R.ipynb` | paralog expression heatmap across 29 samples / 5 organs (**Fig 7D**) |

## Inputs

- **Segdup BEDPE**: `.../biser/hifiasm-041425/segdup_output.bedpe` (step 3) and
  `.../segdup_output_duplicateLinkRemoved.bedpe` (steps 4–5); see `BISER_BEDPE`
  in `config.sh`.
- **Gene annotation** (the new merged annotation, AGAT-normalized):
  `code/command-line-script/annotation-merging/output/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.agat.gff3`
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
