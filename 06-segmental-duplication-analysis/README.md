# 06-segmental-duplication-analysis

Detect segmental duplications (SDs) with BISER across the degu assembly and
related-species genomes.

## Order of operation

| Step | Script | What it does |
|---|---|---|
| 1 | `01_biser_run.sh` | run BISER on one genome → `segdup_output_<species>.bedpe` |
| 2 | `02_biser_submit_batch.sh` | submit step 1 for every genome in `BISER_GENOME_LIST` |

`02` is optional (only needed for the multi-species comparison); run step `01`
directly for a single assembly.

## Downstream processing

BISER's raw `.bedpe` output is refined and analysed elsewhere (not part of this
run step):

- `figure/segdup-investigation/understanding_segdup_output.ipynb` — inspects the
  raw output and trims it to `segdup_output_mod.bedpe`
- `figure/segdup-investigation/regenerate_segdup_overlap.py` /
  `calculate_segdup_geneOverlap_newHiFiAnnotation.ipynb` — SD × gene overlap,
  hotspots, and circos plots

The degu result consumed by the rest of the pipeline is
`.../biser/hifiasm-041425/segdup_output_mod.bedpe` (see `BISER_BEDPE` in
`config.sh`).

## Configuration

Set in `config.sh`: `BISER_BIN`, `BISER_OUT_DIR`, `BISER_THREADS`, `BISER_GC_HEAP`,
`BISER_MAX_EDIT_ERROR`, `BISER_MAX_ERROR`, `BISER_GENOME_DIR`, `BISER_GENOME_LIST`.
The `genome-annotation` conda env (`ENV_GENOME_ANNOTATION`) supplies `samtools`
for indexing.

## Notes

- BISER (v1.4) is a user-level Python install: entry point `~/.local/bin/biser`
  (`BISER_BIN`), package in `~/.local/lib/python3.9/site-packages/biser`.
- BISER requires a **soft-masked** genome; mask repeats before running.
- Scripts are Slurm batch scripts for TSCC/UCSD; adjust the `#SBATCH` lines for
  your cluster.
