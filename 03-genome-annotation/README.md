# 03-genome-annotation

Repeat masking and gene annotation of the degu (*Octodon degus*) assembly.
Two parallel tracks — a **transfer-based** (Liftoff) and a **de-novo**
(Braker3) annotation — are built and then merged into the final annotation.

## Pipeline steps (in order)

| # | Step | Script(s) |
|---|------|-----------|
| 1 | Build de-novo repeat library | `01_repeatModeler.sh` |
| 2 | Repeat masking | `02_repeatMasker.sh` |
| 3 | Transfer-based annotation (Liftoff) | `03-transferBased-geneAnnotation/lift_off.sh` |
| 4 | De-novo annotation (Braker3) | `04-denovoBased-geneAnnotation/` |
| 5 | Merge the two annotations | `05-annotation-merging/` |
| 6 | 3′ UTR annotation + finalize (peaks2utr + AGAT) | `06-annotation-enhancement-peaks2utr-agat/` |

### De-novo track (`04-denovoBased-geneAnnotation/`)

1. `01_submitMultiple_align_transcriptome.sh` — build the HISAT2 index and submit
   one alignment job per RNA-seq SRA accession.
2. `01_align_transcriptome.sh` — align RNA-seq reads + coordinate-sort (array job).
3. `02_submit_brakerScript.sh` — run Braker3 (Singularity) with the aligned BAMs
   + vertebrate protein evidence.
4. `03_blastp_brakerAnnotation.sh` — functional annotation via BLASTP
   (benchmark on a 100-entry subset; the full run used a parallelized script).
5. `04_interproscan_braker.sh` — InterProScan protein-domain annotation (array).
6. `05_annotateGeneLoci_AGAT.txt` — note: annotate gene loci with the
   BLAST/InterProScan results (AGAT).

## Configuration

All paths and settings live in a single `config.sh` at the repository root.
Each script finds and sources it automatically — edit `config.sh` (not the
scripts) to point at your own data, tools, conda environments, and cluster
settings.

## Requirements

- Conda environments: `toolshed-repeatmodeler`, `toolshed-liftoff`,
  `genome-annotation`, `toolshed-peaks2utr`, `toolshed-agat`,
  `toolshed-deeptools`.
- Singularity (`singularitypro/3.11`) for the Braker3 image `braker3.sif`.
- Tools (via `TOOLSHED_DIR`): RepeatModeler/Masker, Liftoff, HISAT2, samtools,
  Braker3, BLAST+, InterProScan.

## Notes

- Scripts are Slurm batch scripts written for TSCC at UCSD; adjust the `#SBATCH`
  scheduler lines for your own cluster.
- `05-annotation-merging/` is a standalone Python tool (see its own README).
- `06-annotation-enhancement-peaks2utr-agat/` reuses the `src/` GFF3 helpers
  from `05-annotation-merging/` (see its own README).
