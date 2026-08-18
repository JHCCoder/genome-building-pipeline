# 02-genome-polishing

Polish the mito + contam-filtered assembly with
[NextPolish2](https://github.com/Nextomics/NextPolish) — a HiFi-aware polishing
pipeline that combines a repetitive-k-mer-masked read mapping (`meryl` +
`winnowmap`) with short-read k-mer datasets (`yak`).

> **Note:** this stage was exploratory. The assembly submitted in the paper
> **skipped polishing** — contamination removal was the only cleaning step
> applied to the final assembly. See `NextPolisher2.doc` for the tool's own
> usage notes.

## Command pipeline (step1–3)

| Step | Script | What it does |
|---|---|---|
| 1 | `01_align_hifi_reads.sh` | `meryl` k-mer DB → `repetitive_k15.txt` → `winnowmap` HiFi mapping → `hifi.map.sort.bam` |
| 2 | `02_kmer_dataset.sh` | `yak` 21-mer and 31-mer datasets (`male403_k21.yak`, `male403_k31.yak`) |
| 3 | `03_nextPolish2_run.sh` | index the BAM, run `nextPolish2` → `genome_chrom_contamRemoved_polished.fa` |

Steps 1 and 2 are independent; step 3 consumes the outputs of both.

## Inputs & outputs

All paths are defined in `config.sh` at the repo root.

| Variable | Meaning |
|---|---|
| `POLISH_ASM` | assembly to polish (mito + contam filtered, pre-scaffolding) |
| `POLISH_HIFI_READS` | HiFi reads for the winnowmap mapping |
| `POLISH_ILLUMINA_R1` / `_R2` | Illumina short reads for the yak k-mer datasets |
| `POLISH_WORK_DIR` | working directory; all intermediates and the polished FASTA land here |

## Notes

- Tools (`meryl`, `winnowmap`, `yak`, `nextPolish2`, `samtools`) are provided by
  the `genome-polishing` conda environment (`ENV_GENOME_POLISHING`).
- Observed runtimes on TSCC (from the original runs, for budgeting): step 1 ≈ 23 h
  (16 threads, ~27 GB), step 2 ≈ 30 min, step 3 ≈ 50 min (8 threads, ~60 GB).
