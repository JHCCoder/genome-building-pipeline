# 06-tandem-repeat-analysis — processing

Tools for detecting tandem repeats / centromeric higher-order repeats, numbered
in the order they are run. This is the **processing** batch (the
analysis/visualization notebooks are added separately).

## Tools (order of operation)

| # | Directory | Tool | What it does |
|---|---|---|---|
| 1 | `01-trf/` | TRF (Tandem Repeats Finder) | Detect tandem repeats: `trf_run.sh` (single), `trf_run_parallelized.sh` + `submit_parallel_trf_run.sh` (split → array), then downstream `run_cd-hit-est.sh`, `blastn_tandem_repeats.sh`, `clustalo_alignment.sh` |
| 2 | `02-centroanno/` | centroAnno | Annotate centromeres, producing monomer templates used by HiCAT |
| 3 | `03-hicat/` | HiCAT | Annotate higher-order repeats (HORs) using the centroAnno templates: `run_HiCAT_chr.sh` (per-chromosome runner), `submit_HiCAT_all.sh` (all 30 chromosomes), `split_chromosome.py` / `split_chr4.py` (chunking), `merge_chunks.py` / `merge_chr25.py` (merge), `resplit_*` (rescue failed chunks) |

## Workflow

1. **TRF** (`01-trf/`) — detect tandem repeats genome-wide; cluster/compare with
   `cd-hit-est` / `blastn` / `clustalo`.
2. **centroAnno** (`02-centroanno/`) — generate centromere monomer templates.
3. **HiCAT** (`03-hicat/`) — annotate centromeric HORs using those templates
   (per-chromosome, chunked, then merged).

## Notes

- Tool binaries and paths are in `config.sh` (repo root): `TRF_BIN`,
  `CENTROANNO_DIR`, `HICAT_BIN`, plus the assembly input paths.
- The HiCAT runner scripts (`run_HiCAT_*`) and split/merge helpers reference
  degu-specific chromosome paths (chr4/chr25 rescues); the generic ones are
  `run_HiCAT_chr.sh` + `submit_HiCAT_all.sh` + `split_chromosome.py` +
  `merge_chunks.py`.
