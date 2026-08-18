# 01-genomic-features (pyCirclize)

Genome-wide circos overview of the degu assembly (chromosomes in ring order,
largest → smallest), built with [pyCirclize](https://github.com/moshi4/pyCirclize).
This is the **no-methylation** version — the CpG-methylation track is removed and
the A/B-compartment and GC-content tracks are re-spaced outward to fill the freed
band.

## Files

| File | What it does |
|---|---|
| `pyCircularize_features_noMethylation.ipynb` | the figure notebook (data prep → plot → savefig) |
| `run_noMethylation_nb.py` | executes the notebook's cells top-to-bottom with the current interpreter |
| `build_noMethylation_nb.py` | one-time derivation script that built the no-methylation notebook from the full `pyCircularize_features.ipynb` (see Notes) |

The notebook renders tracks (inner → outer): A/B compartment, GC content, tRNA,
lncRNA, snRNA, reverse CDS, forward CDS, and contig2scaffold junctions, with
chromosome ticks/labels on the outer ring. Outputs are written to
`degus_genome_circos_genetic_overview_noMethylation.{png,svg,pdf}`.

## Inputs (edit via the `PROJ_ROOT` cell)

| Data | `PROJ_ROOT`-relative path |
|---|---|
| Chromosome sizes (`.fai`) | `data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.hardMasked.mito.fasta.fai` |
| Gene annotation (agat GFF) | `code/command-line-script/annotation-merging/output/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.agat.gff3` |
| GC content (1 kb windows) | `figure/circos-plot/feature-overview/assembly_final.sorted.headerRenamed.chrAssigned.mito.gc_content.bed` |
| A/B compartment scores | `figure/hic-plot/eigenvector_track_dropNA.bed` |
| Contig→scaffold junctions | `figure/circos-plot/feature-overview/agp_final_contig2scaffold.bed` |

## Requirements

- Conda env `python-visualizations` (pyCirclize, matplotlib, pandas, numpy).
- `run_noMethylation_nb.py` needs `nbformat`.

## Notes

- The notebook uses the **correct merged annotation**
  (`…_peaks2utr_sorted.agat.gff3`), substituted in place of the earlier
  `…_peaks2utr_sorted.gff3`. Its `CDS`/`snRNA`/`lnc_RNA`/`tRNA` feature types are
  unchanged, so the gene tracks are unaffected.
- Paths are expressed as `f"{PROJ_ROOT}/…"`; set `PROJ_ROOT` in the first code cell.
- `build_noMethylation_nb.py` reads the full source notebook
  `pyCircularize_features.ipynb` (kept in `figure/circos-plot/feature-overview/`,
  not committed here) to strip the methylation track; it is a derivation record,
  not part of reproducing the figure.
- Generated images (`*.png`/`*.svg`/`*.pdf`) are not committed.
