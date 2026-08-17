# genome-building-pipeline

Code for a de novo genome-building pipeline — assembly & scaffolding, repeat
masking, and genome annotation — for the degu (*Octodon degus*).

## Repository layout

- `01-genome-assembly/` — assembly (hifiasm), Hi-C preprocessing, scaffolding (HapHiC), mitochondrial filtering, sex-chromosome identification
- `02-initial-genome-assembly-evaluation/` — initial assembly QC
- `03-genome-annotation/` — repeat masking, transfer- and de novo-based annotation
- `04-indepth-genome-assembly-assessment/` — in-depth assembly assessment
- `05-tandem-repeat-analysis/` — tandem-repeat analysis
- `06-segmental-duplication-analysis/` — segmental-duplication analysis
- `useful-scripts/` — helper utilities

## Configuration

Edit `config.sh` (repository root) to set paths, conda environments, and cluster
settings for your environment. Every pipeline script finds and sources it
automatically, so you do not need to edit the scripts themselves.

## Notes

- Scripts are Slurm batch scripts written for the TSCC cluster at UCSD; adjust
  the `#SBATCH` scheduler lines for your own cluster.
- Folder names use hyphens (`-`), file names use underscores (`_`).
