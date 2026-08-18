# 07-paralog-loci-read-depth

**Supplemental Figure S7** — HiFi read-coverage comparison between single-copy
genes and paralogs, inside and outside segmental duplications (segdups).

Each locus = gene body **± 5 kb** flank. Five groups are compared:

| Group | Tier | Definition |
|-------|------|-----------|
| BUSCO | BUSCO | single-copy orthologs (glires_odb10) — independent control |
| G1 | T4 | non-paralog gene, **not** overlapping a segdup |
| G2 | T1 | paralog copy (`-lN`/`-dlN`/`-rlN`), **not** overlapping a segdup |
| G3 | T3 | non-paralog gene, overlapping a segdup |
| G4 | T2 | paralog copy (`-lN`/`-dlN`/`-rlN`), overlapping a segdup |

**"Paralog" = paralog copy** — the strict `-lN`/`-dlN`/`-rlN` suffix rule of
`build_paralog_families.py` (04-…/03-paralog-assessment-by-read-depth), so the
G2 + G4 copy counts (2,181) match that analysis. The parent gene is **not**
counted as a paralog; it falls into G1/G3 (non-paralog). **All chromosomes**
(autosomes + chrX/chrY + unplaced scaffolds) are kept, matching that analysis.

The figure has two panels:
- **Panel A** — violin plot of mean HiFi coverage by group (with overlaid boxplot,
  median annotations, and Wilcoxon rank-sum significance brackets).
- **Panel B** — density (KDE) plot of mean HiFi coverage by group.

## Order of operation

| Step | File | What it does |
|------|------|-------------|
| 1 | `01_build_regions.py` | classify genes into tiers; build ±5 kb regions + BUSCO beds |
| 2 | `02_query_coverage.py` | per-base HiFi coverage per region → `region_coverage.tsv` |
| 3 | `03_run_coverage.sh` | `sbatch` wrapper for steps 1–2 (one Slurm job) |
| 4 | `04_coverage_analysis.ipynb` | Panel A violin + Panel B density figure |

## Inputs

- **Gene annotation** (the *new* merged annotation, AGAT-normalized):
  `code/command-line-script/annotation-merging/output/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.agat.gff3`.
  Step 1 reads it indirectly via the gene × segdup overlap table produced from it:
  `07-segmental-duplication-analysis/04_segdup_gene_overlap.ipynb` →
  `hifiasm_gene_segDup_overlapInfo_081726.tsv`. **Re-run that notebook after a
  new annotation before running step 1.**
- **Segdups** (BISER BEDPE, `segdup_output_duplicateLinkRemoved.bedpe`) — used by
  the step-4 overlap notebook, not directly here.
- **HiFi coverage** — per-base coverage track
  (`…_long_read_coverage_basepair_level.per-base.bed.gz`, primary reads mapped to
  the final assembly). Assembly/read-level, not annotation-dependent.
- **BUSCO** — `output/outputs-from-busco-ortholog-alignment/…/run_glires_odb10/full_table.tsv`.
  Assembly-level, not annotation-dependent.
- **Chromosome lengths** — `code/command-line-script/contig-coverage/…_genome_length.txt`.

## Colors

The violin/density groups use the Okabe-Ito colorblind-safe palette (from the
archived notebook `figure/disentangle-haploDuplication-paralogs/_archive/coverage_analysis.orig.ipynb`):

| Group | Tier | Color |
|-------|------|-------|
| G2 | T1 | `#E69F00` |
| G4 | T2 | `#D55E00` |
| G3 | T3 | `#0072B2` |
| G1 | T4 | `#009E73` |
| BUSCO | BUSCO | `#CC79A7` |

## Notes

- Step 2 needs htslib's `tabix` (see the `HTSLIB_MODULES` line in
  `02_query_coverage.py`) and takes ~1–2 h — run it as the Slurm job in step 3,
  not on a login node.
- `region_coverage.tsv` and the `.bed`/`.tsv` intermediate files are generated
  outputs; regenerate them in the numbered order after updating inputs.
- `#SBATCH` scripts target TSCC/UCSD; adjust for your cluster.
