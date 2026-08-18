# 05-structural-error-detection-Inspector

Detection of assembly structural errors and small-scale errors with
[Inspector](https://github.com/ChongLab/Inspector), a reference-free long-read
assembly evaluator.

## Files

| File | What it does |
|---|---|
| `01_inspector_run.sh` | run Inspector on the scaffolded assembly with HiFi reads (reference-free; reference-based mode commented) |
| `METHODS.txt` | methods-section description (commands, inputs, metrics, citation) |

## Outputs (in `$INSPECTOR_OUT_DIR`)

Inspector writes, per run:
- `summary_statistics` — contig N50, mapping rate, depth, structural errors
  (expansion/collapse/haplotype-switch/inversion), small-scale errors, QV.
- `structural_error.bed` / `small_scale_error.bed` — the error calls.
- `Inspector.log`, `valid_contig.fa`, read-to-contig BAM, etc.

The report files (`summary_statistics`, `structural_error*.bed`,
`small_scale_error*.bed`, `Inspector.log`, `contig_length_info`) are copied to
`100-outputs-genome-information/Inspector/<run>/`; the large BAM/FASTA files are
not stored in the repo.

## Configuration

Inputs and the conda env live in `config.sh` (`ENV_INSPECTOR`, `INSPECTOR_BIN`,
`INSPECTOR_ASM`, `INSPECTOR_HIFI_READS`, `INSPECTOR_REF`, `INSPECTOR_OUT_DIR`).

## Notes

- `#SBATCH` targets TSCC/UCSD; adjust for your cluster.
- The `-r` argument is unquoted so the wildcard in `INSPECTOR_HIFI_READS` expands.
