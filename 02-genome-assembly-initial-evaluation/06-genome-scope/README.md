# 06-genome-scope

Estimate genome size, heterozygosity, and repeat content from the Illumina WGS
k-mer spectrum with [GenomeScope](http://genomescope.org/).

`genome_scope.sh` runs KMC to count k-mers and write a frequency histogram
(`*.histo`). **The histogram is not fitted here** — instead, upload it to the
GenomeScope web interface to get the estimate.

## Workflow

1. `genome_scope.sh` counts k-mers with KMC (`-k 21`) from the reads in
   `GENOMESCOPE_READS` and writes `<out>/reads<cs>.histo`.
2. Upload that `.histo` file to the GenomeScope web page
   ([genomescope.org](http://genomescope.org/), or
   [GenomeScope 2.0](http://qb.cshl.edu/genomescope/) for the newer model),
   set the matching `k`, and read off the estimated genome size,
   heterozygosity, and repeat fraction.

## Run parameters (top of `genome_scope.sh`)

| Variable | Meaning |
|---|---|
| `KMER` | k-mer size (21) |
| `MIN_COVERAGE` / `MAX_COVERAGE` | KMC `-ci` / `-cs` coverage cutoffs |
| `MEM` | KMC max memory (GB) |
| `THREADS` | KMC thread count |
| `SAMPLE` | label used in the output directory name |

## Configuration

Paths and KMC binaries are set in `config.sh`: `KMC_BIN`, `KMC_TOOLS_BIN`,
`GENOMESCOPE_READS`, `GENOMESCOPE_OUT_DIR`, `GENOMESCOPE_TEMP_DIR`.

## Notes

- Scripts are Slurm batch scripts for TSCC/UCSD; adjust the `#SBATCH` lines for
  your cluster.
