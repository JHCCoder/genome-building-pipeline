#!/bin/bash
#SBATCH -J polish_align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 2-00:00:00
#SBATCH --mem=80G
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

mkdir -p "$POLISH_WORK_DIR"
cd "$POLISH_WORK_DIR"

# Build a k-mer frequency DB of the assembly and extract the repetitive k-mers (winnowmap -W)
meryl count k=15 "$POLISH_ASM" output hifi0414_aggressivePurge_mitoContamRemoved.merylDB
meryl print greater-than distinct=0.9998 hifi0414_aggressivePurge_mitoContamRemoved.merylDB > repetitive_k15.txt

# Map HiFi reads to the assembly (winnowmap map-pb), sort to BAM
winnowmap -t 16 -W repetitive_k15.txt -ax map-pb "$POLISH_ASM" $POLISH_HIFI_READS | samtools sort -o hifi.map.sort.bam -
