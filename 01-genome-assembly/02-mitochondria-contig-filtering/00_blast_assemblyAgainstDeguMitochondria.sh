#!/bin/bash
#SBATCH -J mito_chrom_blast
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 03:00:00
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

assembly_alias="hifiasm_041425"
assembly="$MITO_BLAST_ASSEMBLY"
mito_chrom="$MITO_FASTA"

makeblastdb -in "$assembly" -dbtype nucl
blastn -query "$mito_chrom" -db "$assembly" \
  -outfmt "6 qseqid sseqid pident qcovs length gapopen evalue bitscore mismatch" \
  -evalue 1e-5 > "${assembly_alias}_mito_blast.out"
