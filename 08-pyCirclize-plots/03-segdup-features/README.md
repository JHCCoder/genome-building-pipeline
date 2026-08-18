# 03-segdup-features (pyCirclize)

Genome-wide circos of the degu assembly showing segmental-duplication (segdup)
features, built with [pyCirclize](https://github.com/moshi4/pyCirclize). Renders
(inner → outer): segdup density at gene-level (blue), segdup density at bp-level
(gold), segdup links (BISER BEDPE, subsampled), contig→scaffold junctions
(dotted), and chromosome ticks/labels on the outer ring. Outputs to
`degus_genome_circos_segdup_overview.{png,svg,pdf}`.

## Files

| File | What it does |
|---|---|
| `pyCircularize_segdup_overview.ipynb` | the figure notebook (imports → data prep → plot → savefig) |

## Inputs (edit via the `PROJ_ROOT` cell)

| Data | `PROJ_ROOT`-relative path |
|---|---|
| Chromosome sizes (`.fai`) | `data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.hardMasked.mito.fasta.fai` |
| Segdup density, gene-level | `figure/segdup-investigation/segdup_density_rolling_5000000_100000_geneLevel.tsv` |
| Segdup density, bp-level | `figure/segdup-investigation/segdup_density_rolling_5000000_100000_bpLevel.tsv` |
| Segdup links (BISER BEDPE) | `code/command-line-script/genome-annotation/biser/hifiasm-041425/segdup_output_mod.bedpe` (== `BISER_BEDPE` in `config.sh`) |
| Contig→scaffold junctions | `figure/circos-plot/feature-overview/agp_final_contig2scaffold.bed` |

## Requirements

- Conda env `python-visualizations` (pyCirclize, matplotlib, pandas, numpy).

## Notes

- The source notebook (`figure/circos-plot/segdup-feature/pyCircularize.ipynb`,
  148 cells) was reduced to the 12 cells that produce this figure; tutorial
  ("Load Genbank …", `load_*_example_dataset`), GPT/exploratory, and duplicate
  cells were removed, and the two `from pycirclize import Circos` /
  `to_rgba` duplicates were deduplicated.
- The segdup links are **subsampled to 10%** (`bedpe.sample(frac=0.1, random_state=28)`)
  for plotting speed; a `_full` (100%) variant was also produced in the source but
  is not reproduced here.
- **Stale input caveat:** `…_geneLevel.tsv` (Oct 2025) was computed from the
  pre-merge annotation; the gene×segdup overlap was since regenerated from the
  merged AGAT GFF (`hifiasm_gene_segDup_overlapInfo_092525.tsv`, Aug 2026). If
  gene-overlap changed materially, regenerate the gene-level density before reuse.
  The bp-level density and BEDPE are annotation-independent.
- The source notebook also carried six assembly/reference variants
  (`pyCircularize_{verkko,purged,mouse,human,guineaPig}.ipynb`) and a separate
  "annotation-change" figure (nc/gr/ga gene tracks → `degus_genome_circos_annotation_change_overview.png`);
  those are out of scope here.
- Generated images (`*.png`/`*.svg`/`*.pdf`) are not committed.
