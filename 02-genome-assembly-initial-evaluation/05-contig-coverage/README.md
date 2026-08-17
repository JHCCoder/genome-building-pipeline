# 05-contig-coverage

Map short (Illumina WGS) and long (HiFi) reads to an assembly and compute
read-depth coverage in fixed windows, per-contig, and genome-wide
(`bedtools coverage`).

One script (`map_and_coverage.sh`) covers both read types. Edit the run
parameters at the top of the script; all paths and conda environments live in
`config.sh` at the repo root.

## Modes (two independent toggles)

| `MULTIMAP` | `FILTER` | `samtools view` flags | Meaning |
|---|---|---|---|
| `yes` | `no` | (none) | keep every alignment (multi-mapping, unfiltered) |
| `no`  | `no` | `-F 0x100` | unique — drop secondary alignments |
| `yes` | `yes` | `-F 0x800 -q MAPQ` | multi-mapping, filtered (drop chimeric + low-MAPQ) |
| `no`  | `yes` | `-F 0x900 -q MAPQ` | unique + filtered (strictest) |

`READ_TYPE` picks the aligner:

- `short` → `bwa-mem2 mem` (multi-mapping adds `-a`, since bwa does not emit
  secondary alignments by default)
- `long` → `minimap2 -ax map-hifi` (unique adds `--secondary=no`)

## Run parameters (top of `map_and_coverage.sh`)

| Variable | Meaning |
|---|---|
| `READ_TYPE` | `short` or `long` |
| `MULTIMAP` | `yes` (multi-mapping) / `no` (unique) |
| `FILTER` | `yes` (drop supplementary + MAPQ) / `no` (no filtering) |
| `MAPQ` | mapping-quality threshold (only used when `FILTER=yes`) |
| `WINDOW_SIZE` | coverage window size in bp (default 15000) |
| `ASSEMBLY_ALIAS` | label used in output filenames |

## Outputs (written to `COVERAGE_OUT_DIR`)

- `<alias>_<read>_read_<mode>_alignment_sorted.bam` (+ `.bai`) — mapped BAM
- `coverage_<read>_read_<alias>_<mode>_<N>bp_windows.tsv` — mean coverage per window
- `coverage_<read>_read_<alias>_<mode>_contig.tsv` — mean coverage per contig
- `coverage_<read>_read_<alias>_<mode>_whole_assembly.tsv` — genome-wide mean

## Configuration

Paths and conda environments are set in `config.sh`:
`COVERAGE_ASSEMBLY`, `COVERAGE_SHORT_R1`/`_R2`, `COVERAGE_HIFI_READS`,
`COVERAGE_OUT_DIR`, `ENV_BULK_HIC` (bwa-mem2 + samtools), `ENV_ASSESSMENT`
(minimap2 + samtools), `ENV_BEDTOOLS` (bedtools).

## Notes

- The BWA-MEM2 index and samtools `.fai` are written next to the assembly
  (guarded by existence checks).
- Scripts are Slurm batch scripts for TSCC/UCSD; adjust the `#SBATCH` lines for
  your cluster.
