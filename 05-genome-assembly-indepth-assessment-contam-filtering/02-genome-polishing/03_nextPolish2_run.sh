#!/bin/bash
#SBATCH -J polish_run
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 2-00:00:00
#SBATCH --mem=200G
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
conda activate "$ENV_GENOME_POLISHING"

cd "$POLISH_WORK_DIR"

# Index the HiFi BAM (from step1), then polish with NextPolish2 using the k-mer datasets (step2)
samtools index -@ 8 hifi.map.sort.bam
nextPolish2 -t 8 hifi.map.sort.bam "$POLISH_ASM" male403_k21.yak male403_k31.yak > genome_chrom_contamRemoved_polished.fa
