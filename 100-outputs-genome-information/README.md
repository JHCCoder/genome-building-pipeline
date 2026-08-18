# 100-outputs-genome-information

Actual output / result files for the de-novo degu (*Octodon degus*) genome
assembly (OctDeg2.0), kept separate from the pipeline code (stages `01–07`).
Only reasonably-sized reports and text files are stored here; very large files
are deposited in a data repository and referenced by accession.

## The genome assembly itself is too large to commit

The final assembly FASTA
(`assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta`, **3.4 Gb**,
288 scaffolds, N50 ~125 Mb) is **not stored in this repository** — it is far
too large for git. It can be shared **by email on request** and will be
deposited to **NCBI GenBank** when the manuscript is published (accession to be
added to the Data Availability section). Everything else small enough to store
is deposited here.

## Quick QC summary

| Metric | Tool | Value |
|---|---|---|
| Gene completeness | BUSCO v5.7.1 (glires_odb10, n=13,798) | **C:98.7%** [S:90.1%, D:8.6%], F:0.7%, M:0.6% |
| K-mer completeness | Merqury | **87.72%** |
| Consensus quality (QV) | Merqury | **QV 32.14** |
| Assembly size / contiguity | BUSCO assembly stats | 3.40 Gb, 288 scaffolds, N50 125 Mb |
| Structural errors | Inspector v1.3.1 (HiFi) | **62** (44 expansion, 11 collapse, 3 haplotype switch, 4 inversion) |
| Small-scale error rate | Inspector v1.3.1 (HiFi) | **17.18 / Mb** |
| Consensus quality (QV, HiFi) | Inspector v1.3.1 | **QV 45.96** |

## Contents

| Directory | Analysis | Key files |
|---|---|---|
| `annotation/` | final merged gene annotation (Liftoff + Braker3/AUGUSTUS) | `hifiasm-041425-denovoEnhanced_peaks2utr_sorted.agat.gff3.gz` |
| `segdup-biser/` | segmental duplications | BISER `segdup_output_duplicateLinkRemoved.bedpe.gz` (primary call set), `segdup_output_mod.bedpe.gz`, `…_elem.txt.gz` |
| `centroAnno/` | centromere annotation (per chromosome + genome-wide) | centroAnno `chr*_HORs.fa`, `chr*_monomerTemplates.fa`, `chr*_{hor,}decomposedResult.csv.gz`, plus `cautils-*/` genome-wide `HORs.bed`, `repeat_regions.bed`, `top10_repeats.bed` |
| `hicat/` | higher-order-repeat arrays (chr4 & chr25) | HiCAT `hicat_chr{4,25}_{density,hor_arrays,types}.bed` |
| `trf/` | tandem-repeat finder | `repeat_df_degu.tsv.gz` (full degu TRF output, 350 MB → 46 MB), `tanRepeatDensity.tsv` / `tanRepeatDensity.merged.tsv` (1 Mb-window counts / rolling density), `repeatDensity…rollingWindow.tsv` (all repeats), plus 3 locus/family `.dat` runs |
| `fcs-gx/` | cross-species contamination screen | `fcs_gx.contam.action_report.txt` |
| `fcs-adaptor/` | adaptor contamination screen | `adaptor_report.txt` |
| `busco/` | gene-content completeness | `short_summary.specific.glires_odb10…txt` (+ `.json`) |
| `merqury/` | k-mer completeness + consensus QV | `*.completeness.stats`, `*.qv` |
| `Inspector/` | assembly structural-error detection | Inspector v1.3.1 reports |

## Notes

- **TRF** — the full degu TRF output is `trf/repeat_df_degu.tsv.gz` (parsed from
  the whole-genome `.dat`, one row per tandem repeat; 350 MB → 46 MB gzipped).
  The whole-genome `.dat` (3.8 Gb) and the per-species parsed outputs
  (`repeat_df_<species>.tsv`, 121–573 Mb each) are too large to commit and live
  in `code/command-line-script/genome-annotation/trf-tandem-repeat/`. Summaries
  are in `trf/` (`tanRepeatDensity.tsv` = TRF count per 1 Mb window,
  `tanRepeatDensity.merged.tsv` = rolling density) and in the centroAnno
  TRF-derived units (`centroAnno/` `*_monomerTemplates.fa` +
  `cautils-*/repeat_regions.bed` + `cautils-*/top10_repeats.bed`). `trf/` also
  holds 3 standalone locus/family `.dat` runs.
- Large text files are stored gzipped (`*.gz`) — `gunzip` to decompress.
- `segdup-biser/segdup_output_duplicateLinkRemoved.bedpe.gz` is the primary
  (reciprocal-duplicate-removed) segmental-duplication call set.
- `hicat/` covers the **two chromosomes with the largest HOR arrays — chr4 and
  chr25**.
