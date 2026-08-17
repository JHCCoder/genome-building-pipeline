# 02-initial-genome-assembly-evaluation

Initial quality assessment of the finished assembly: Merqury (k-mer-based
completeness/consensus), BUSCO (gene-content completeness), and basic FASTA
metrics (N50/L50, GC, contig-size distribution).

## Pipeline steps (in order)

| # | Step | Script(s) |
|---|------|-----------|
| 1 | Build the Illumina WGS k-mer database (Merqury prerequisite) | `01_meryl_db.sh` |
| 2 | Run Merqury against the assembly | `02_merqury.sh` |
| 3 | Run BUSCO (genome mode, `glires_odb10` lineage) | `03_busco_ortholog_alignment.sh` |
| 4 | Compute basic FASTA metrics (N50/L50/GC) | `04_fasta_metrics.py` |

BUSCO completeness/duplication metrics (C/S/D/F/M, n) are read directly from the
resulting `short_summary.specific.<lineage>.<prefix>.txt` file; no BUSCO plotting
is performed.

## Configuration

All paths and settings live in a single `config.sh` at the repository root. Each
script finds and sources it automatically — edit `config.sh` (not the scripts) to
point at your own data, tools, conda environments, and cluster settings.

## Requirements

- Conda environments: `toolshed-merqury` (Merqury + meryl), `toolshed-busco` (BUSCO).
- `04_fasta_metrics.py` needs Python 3 with Biopython.

## Notes

- Scripts are Slurm batch scripts written for TSCC at UCSD; adjust the `#SBATCH`
  scheduler lines for your own cluster.
