# 08-genome-comparison-with-octDeg1

**Figure 3** — comparison of the two genome assemblies (OctDeg1.0 reference vs
the new OctDeg2.0 / hifiasm-041425), showing that the new assembly adds
repetitive regions and captures more epigenetic / single-cell signal.

Four panels (revised for Genome Research, October 2026; the earlier seven-panel
layout merged B+C, D+E and F+G), from four data sources:

| Panel | Description | Data source | Script | Size |
|-------|-------------|-------------|--------|------|
| A | Whole-genome alignment dotplot (minimap2 `-x asm5`) | `.paf` | `01_alignment_scatterplot.py` | 110 × 95 mm |
| B | Repetitive / non-repetitive OctDeg2.0 sequence (Gb), aligned vs not aligned to OctDeg1.0 | aligned BED × RepeatMasker | `02_repeat_overlap_barplots.py` | 58 × 95 mm |
| C | Bulk Hi-C: cis contacts (top) and cis/trans ratio (bottom), 4 samples | `cis_trans_bhic_metrics.txt` | `03_bhic_cistrans_barplots.py` | 84 × 60 mm |
| D | scMultiome: feature linkages (top) and median genes per cell (bottom), 4 samples, matched gene models | `feature_linkage.txt`, `gene_per_cell.txt` | `04_scmultiome_metric_barplots.py` | 84 × 60 mm |

Layout on the 175 mm page: A | B in the first row, C | D in the second.
Every script draws at printed size (8 pt Helvetica/Nimbus Sans, no titles) and
writes PDF + SVG + PNG (600 dpi) to `figure3/` (`$FIG3_OUT_DIR`); an existing
output is kept as `*.pre-<timestamp>` before it is replaced.

## Order of operation

| Step | File | What it does |
|------|------|-------------|
| 1 | `01_alignment_scatterplot.py` | parse + filter the PAF, orient contigs, plot Panel A |
| 2 | `02_repeat_overlap_barplots.py` | aligned × repeat basepair overlap (cached in `figure3/figure3B_repeat_alignment_bp.tsv`; `--recompute` to rebuild); plot Panel B |
| 3 | `03_bhic_cistrans_barplots.py` | cis contacts + cis/trans ratio (Panel C) |
| 4 | `04_scmultiome_metric_barplots.py` | feature linkages + median genes per cell (Panel D) |

Run with the `python-visualizations` env, e.g. `python 01_alignment_scatterplot.py`.
Each script accepts `--width-mm` / `--height-mm` if the layout changes.

## Inputs

Large inputs are read through `config.sh` at the repo root (`FIG3_PAF`,
`FIG3_ALIGNED_BED`, `FIG2B_REPEAT_GFF`, `FIG3_CHROM_FAI`); the four small
per-sample summary tables are bundled in `data/`.

- **PAF** (Panel A): `figure/align-octDeg1-hifiasm-minimap2/alignment_octDeg1_hifi041425.paf`
  — minimap2 `-x asm5`, query = OctDeg2.0, target = OctDeg1.0.
- **Aligned BED** (Panel B): `figure/align-octDeg1-hifiasm-minimap2/hifi041425_alignment.bed`
  — query-side aligned intervals (OctDeg2.0 side).
- **RepeatMasker GFF** (Panel B): `figure/circos-plot/feature-overview/assembly_final.sorted.headerRenamed.fasta.out.chr.gff`
  — repetitive regions of the final OctDeg2.0 assembly.
- **Chromosome sizes** (Panel B): `…/assembly_final.sorted.headerRenamed.chrAssigned.fasta.fai`.
- **bHiC metrics** (Panel C): `data/cis_trans_bhic_metrics.txt`
  (`Sample  Cis_Trans  Cis_Contact  Read_mapped  Pairs_mapped  Assembly`).
- **scMultiome metrics** (Panel D): `data/feature_linkage.txt`,
  `data/gene_per_cell.txt` (and `data/fragment_per_cell.txt`, an extra metric).

## Samples

- Bulk Hi-C (Panel C): `degu_060302`, `degu_120106`, `degu_060303`, `degu_040604`
  (shown without the `degu_` prefix).
- scMultiome (Panel D): `181_PFC`, `181_dHIP`, `6997_dHIP`, `7000_PFC`. OctDeg1.0 values come from
  the original OctDeg1 runs; OctDeg2.0 values from the matched-annotation reruns
  (`05_build_liftoff_only_gtf.py` → `06_cellranger_arc_mkref_liftoff.sh` →
  `07_cellranger_arc_count_octDeg2.sh` → `08_update_fig3d_octdeg2_metrics.py`; libraries per
  sample in `data/fig3d_multiome_libraries.tsv`)
  (shown with a space instead of `_`; editable via `TARGET_SAMPLES`).

## Methods notes

- **Panel A** keeps alignments with mapping quality ≥ 60 and alignment length
  ≥ 10 kb, and only those on each OctDeg1.0 contig's dominant strand
  (syntenic). The y-axis orders OctDeg1.0 contigs by their primary OctDeg2.0
  chromosome and orients them by their dominant strand, so collinear alignments
  fall on the diagonal. Current numbers: 202,707 raw alignments → 71,301 after
  the mapq/length filter → 68,607 syntenic, from **674** OctDeg1.0 contigs
  (2.96 Gb); median identity 99.2%. (The old panel title said "1,735 aligned
  contigs"; that number was stale.) Chromosomes are labelled by number under
  the x-axis and shaded in alternating tones instead of a colour key.
- **Panel B** merges overlapping aligned and repetitive intervals and computes
  exact basepair overlap with a two-pointer sweep (pure pandas/numpy). One bar
  per repeat class (height = sequence in Gb); the dark part is aligned. Above
  each bar: % of the class aligned (formerly panel C; 91.9% non-repetitive,
  57.4% repetitive). Inside: the class's share of all aligned sequence
  (formerly panel B; 61.8% / 38.2%). The RepeatMasker annotation is on the
  *final* assembly while the alignment is against hifiasm-041425, so chromosome
  sizes are taken from the final assembly `.fai` for a consistent
  non-repetitive denominator.
- **Panels C/D** are grouped barplots of pre-computed per-sample summaries; no
  alignment is re-run here.

## Colors

Both assemblies use the same two colors throughout, defined once in `config.sh`:
OctDeg1.0 = `#999999` (grey, `FIG3_COLOR_OCTDEG1`), OctDeg2.0 = `#0072B2`
(Okabe–Ito blue, `FIG3_COLOR_OCTDEG2`). Panel A alternates OctDeg2.0 blue and a
light tint by chromosome; panel B uses the same blue (aligned) and light tint
(not aligned). These replace the earlier pale cyan/yellow pair.

## Notes

- Panel D shows the **median** genes per cell (Cell Ranger metric), not the mean.
- Feature linkages and genes per cell are **scMultiome** (10x Multiome)
  metrics, not Hi-C metrics.
- The `figure3*_*.png` files directly in this folder are the pre-revision
  seven-panel outputs; the current panels are in `figure3/`.
