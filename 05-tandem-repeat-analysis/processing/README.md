# 05-tandem-repeat-analysis — processing

Tools for detecting tandem repeats / centromeric higher-order repeats. This is
the **processing** batch (the analysis/visualization notebooks are added
separately).

## Tools

| Directory | Tool | What it does |
|---|---|---|
| `trf/` | TRF (Tandem Repeats Finder) | Detect tandem repeats: `trf_run.sh` (single), `trf_run_parallelized.sh` + `submit_parallel_trf_run.sh` (split → array), plus downstream `run_cd-hit-est.sh`, `blastn_tandem_repeats.sh`, `clustalo_alignment.sh` |
| `centroanno/` | centroAnno | Annotate centromeres (produces monomer templates used by HiCAT) |
| `hicat/` | HiCAT | Higher-order repeat / HOR annotation: `run_HiCAT_chr.sh` (per-chromosome runner), `submit_HiCAT_all.sh` (all 30 chromosomes), `split_chromosome.py` / `split_chr4.py` (chunking), `merge_chunks.py` / `merge_chr25.py` (merge), `resplit_*` (rescue failed chunks) |

## Workflow

1. **centroAnno** — generate centromere monomer templates from the assembly.
2. **HiCAT** — use those templates to annotate higher-order repeats
   (per-chromosome, chunked, then merged).
3. **TRF** — detect tandem repeats; downstream `cd-hit-est`/`blastn`/`clustalo`
   cluster and compare the repeat sequences.

## Notes

- Tool binaries and paths are in `config.sh` (repo root): `TRF_BIN`,
  `CENTROANNO_DIR`, `HICAT_BIN`, plus the `MITO_ASSEMBLY`/`HICAT` input paths.
- The HiCAT runner scripts (`run_HiCAT_*`) and split/merge helpers reference
  degu-specific chromosome paths (chr4/chr25 rescues); the generic ones are
  `run_HiCAT_chr.sh` + `submit_HiCAT_all.sh` + `split_chromosome.py` +
  `merge_chunks.py`.
