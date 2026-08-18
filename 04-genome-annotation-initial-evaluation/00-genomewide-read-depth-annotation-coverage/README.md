# 00-genomewide-read-depth-annotation-coverage

Genome-wide per-chromosome "track" figures: read-depth coverage (primary
alignments) and annotation-density, each drawn as one horizontal track per
chromosome. Pulled from `figure/track-plot/`.

## Files

| File | Figure |
|---|---|
| `01_read_depth_tracks_primaryReads.ipynb` | HiFi (long) + Illumina (short) primary-read coverage, unmasked vs hard-masked, 150 kb bins → `coverage_genome_primaryReadOnly.png` and `coverage_genome_rolling_primaryRead.png` |
| `02_annotation_density_tracks.ipynb` | Liftoff + Braker3 gene density and RepeatMasker repeat density, 1 Mb windows → `density_count_annotation.png` and `density_count_annotation_rolling.png` |

## Inputs (edit via the `PROJ_ROOT` cell)

- **Coverage** (`01_…`): 4 `.tsv` files in
  `code/command-line-script/contig-coverage/2-primaryreads-coverage/`
  (long/short × unmasked/hard-masked, 15 kb windows; hard-masked long reads use
  the `…_fromCalculateScript.tsv` variant).
- **Chromosome sizes**: `data/…/assembly_final.sorted.headerRenamed.chrAssigned.hardMasked.fasta.fai`.
- **Annotations** (`02_…`):
  - Liftoff: `output/outputs-from-liftoff/hifiasm-041425-scaffolded-chrAssigned/hifiasm-041425-scaffolded-chrAssigned.gff`
  - Braker3: `code/command-line-script/genome-annotation/annotate-braker3-results/annotate-gff/braker.uniprotBlast.interproscan.gff/braker.gff3`
  - Repeats: `figure/circos-plot/feature-overview/assembly_final.sorted.headerRenamed.fasta.out.chr.gff`

## Requirements

- Python with pandas, matplotlib, numpy, natsort (e.g. conda env `genome-assembly`).

## Notes

- Both notebooks were reduced from their exploratory sources (55 and 33 cells →
  11 and 17 cells): removed single-locus lookups, one-off `head()`/`shape`/
  `display` checks, unused imports (`Line2D`, `figure`, `seaborn`, `mpatches`),
  and the unused `calculate_density_from_gff` helper. The retained plot cells are
  verbatim.
- The annotation notebook plots the two *sources* (Liftoff transfer + Braker3
  de-novo) separately, not the merged AGAT GFF; the Liftoff/Braker3/repeat GFFs
  are the current ones.
- Generated figures (`*.png`) are not committed.
