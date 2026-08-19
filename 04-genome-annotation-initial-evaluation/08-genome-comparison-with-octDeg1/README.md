# 08-genome-comparison-with-octDeg1

**Figure 3** — comparison of the two genome assemblies (OctDeg1.0 reference vs
the new OctDeg2.0 / hifiasm-041425), showing that the new assembly adds
repetitive regions and captures more epigenetic / single-cell signal.

Seven panels, from four data sources:

| Panel | Description | Data source | Script |
|-------|-------------|-------------|--------|
| A | Whole-genome alignment scatterplot (minimap2 `-x asm5`) | `.paf` | `01_alignment_scatterplot.py` |
| B | % of aligned regions within repetitive / non-repetitive regions | aligned BED × RepeatMasker | `02_repeat_overlap_barplots.py` |
| C | % of repetitive / non-repetitive regions aligned to OctDeg1.0 | aligned BED × RepeatMasker | `02_repeat_overlap_barplots.py` |
| D | Cis contact (bulk Hi-C, 4 samples) | `cis_trans_bhic_metrics.txt` | `03_bhic_cistrans_barplots.py` |
| E | Cis/Trans ratio (bulk Hi-C, 4 samples) | `cis_trans_bhic_metrics.txt` | `03_bhic_cistrans_barplots.py` |
| F | Feature linkages (scMultiome, 4 samples) | `feature_linkage.txt` | `04_scmultiome_metric_barplots.py` |
| G | Median genes per cell (scMultiome, 4 samples) | `gene_per_cell.txt` | `04_scmultiome_metric_barplots.py` |

## Order of operation

| Step | File | What it does |
|------|------|-------------|
| 1 | `01_alignment_scatterplot.py` | parse + filter the PAF, orient contigs, plot Panel A |
| 2 | `02_repeat_overlap_barplots.py` | bedtools overlap of aligned × repeat; plot Panels B, C |
| 3 | `03_bhic_cistrans_barplots.py` | barplots of cis contact / cis-trans ratio (Panels D, E) |
| 4 | `04_scmultiome_metric_barplots.py` | barplots of feature linkage / genes per cell (Panels F, G) |

## Inputs

Large inputs are referenced by absolute path (edit `PROJ_ROOT` at the top of
each script); the four small per-sample summary tables are bundled in `data/`.

- **PAF** (Panel A): `figure/align-octDeg1-hifiasm-minimap2/alignment_octDeg1_hifi041425.paf`
  — minimap2 `-x asm5`, query = OctDeg2.0, target = OctDeg1.0.
- **Aligned BED** (Panels B/C): `figure/align-octDeg1-hifiasm-minimap2/hifi041425_alignment.bed`
  — merged query-side aligned intervals (OctDeg2.0 side).
- **RepeatMasker GFF** (Panels B/C): `figure/circos-plot/feature-overview/assembly_final.sorted.headerRenamed.fasta.out.chr.gff`
  — repetitive regions of the final OctDeg2.0 assembly.
- **Chromosome sizes** (Panels B/C): `…/assembly_final.sorted.headerRenamed.chrAssigned.fasta.fai`.
- **bHiC metrics** (Panels D/E): `data/cis_trans_bhic_metrics.txt`
  (`Sample  Cis_Trans  Cis_Contact  Read_mapped  Pairs_mapped  Assembly`).
- **scMultiome metrics** (Panels F/G): `data/feature_linkage.txt`,
  `data/gene_per_cell.txt` (and `data/fragment_per_cell.txt`, an extra metric).

## Samples

- Bulk Hi-C (Panels D/E): `degu_060302`, `degu_060303`, `degu_120106`, `degu_040604`.
- scMultiome (Panels F/G): `181_PFC`, `181_dHIP`, `6997_dHIP`, `7000_PFC`
  (editable via `TARGET_SAMPLES` in `04_scmultiome_metric_barplots.py`).

## Methods notes

- **Panel A** keeps alignments with mapping quality ≥ 60 and alignment length
  ≥ 10 kb, and (optionally) only the syntenic strand per contig. The y-axis
  orders OctDeg1.0 contigs by their primary OctDeg2.0 chromosome and orients
  them by their dominant strand, so collinear alignments fall on the diagonal.
- **Panels B/C** merge overlapping aligned and repetitive intervals, compute
  exact basepair overlap with `bedtools intersect -wo`, then express it against
  two denominators (total aligned bp for B; total repetitive / non-repetitive
  genome bp for C). The RepeatMasker annotation is on the *final* assembly while
  the alignment is against hifiasm-041425, so chromosome sizes are taken from
  the final assembly `.fai` for a consistent non-repetitive denominator.
- **Panels D–G** are grouped barplots of pre-computed per-sample summaries; no
  alignment is re-run here.

## Colors

Both assemblies use the same two colors throughout: `octDeg1` = `#33ffff`
(cyan), `octDeg2` = `#ffcc00` (yellow). Panel A colors chromosomes by the
`turbo` colormap.

## Notes

- The figure caption says "mean genes per cell" for Panel G, but the Cell Ranger
  metric (and `gene_per_cell.txt`) is the **median** genes per cell; the script
  uses the median value and labels the axis accordingly.
- The caption's "from 4 Hi-C samples" for Panels F/G is a copy-editing slip —
  feature linkages and genes per cell are **scMultiome** (10x Multiome) metrics,
  not Hi-C metrics.
- `02_repeat_overlap_barplots.py` needs `pybedtools` + the `bedtools` binary.
  Scripts 3–4 read their inputs from `data/`; scripts 1–2 read large inputs via
  `PROJ_ROOT` absolute paths.
