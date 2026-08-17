#!/bin/bash
#SBATCH -J sry_blast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH --mem=64G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_GENOME_ANNOTATION"

assembly_alias="hifiasm_041425_scaffolded_curated_masked"
assembly="$SRY_BLAST_ASSEMBLY"
sry_gene="$SRY_GENE"

# If a blast db has not been built for this assembly yet, uncomment:
# makeblastdb -in "$assembly" -dbtype nucl

blastn -query "$sry_gene" -db "$assembly" -outfmt 5 -out results_for_kablammo.xml
