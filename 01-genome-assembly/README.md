# 01-genome-assembly

De novo assembly and scaffolding of the degu (*Octodon degus*) genome.

## Pipeline steps (in order)

| # | Step | Script(s) |
|---|------|-----------|
| 0 | Hi-C read preprocessing (trim → align → pairs → contact matrix) | `00-bHICfiles-preprocessing/` (submission + worker) |
| 1 | Assembly with hifiasm (HiFi + Hi-C mode) | `01_genome_assembly_hifiasm_HiCMode.sh` |
| 2 | Mitochondrial contig filtering | `02-mitochondria-contig-filtering/` |
| 3 | Hi-C scaffolding with HapHiC | `03_genome_scaffolding_haphic.sh` |
| 4 | Sex-chromosome identification (synteny / coverage / kmers / SRY) | `04-sexChr-identification/` |
| 5 | Scaffold curation & chromosome naming | `05_curate_fasta.txt` (manual step) |

## Configuration

All paths and settings live in a single `config.sh` at the repository root.
Each script finds and sources it automatically, so you only edit `config.sh`
(not the scripts) to point at your own data, tools, conda environments, and
cluster settings.

## Requirements

- Conda environments: `bulk-HiC-processing`, `haphic`, `genome-annotation`,
  `toolshed-DiscoverY`, `toolshed-jcvi`.
- Tools (via conda or `TOOLSHED_DIR`): hifiasm, HapHiC, bwa/bwa-mem2, samtools,
  pairtools, cooler, trim_galore, blast+, gffread, jcvi, DiscoverY.

## Notes

- Scripts are Slurm batch scripts written for TSCC at UCSD. Adjust the
  `#SBATCH` scheduler lines for your own cluster.
- `04-sexChr-identification/01-synteny-to-sex-chr/` (JCVI/MCscan synteny to sex
  chromosomes) includes its helper scripts; it also needs the mouse GRCm39
  reference genome/GFF downloaded from NCBI (see its script header).
