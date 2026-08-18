#!/bin/bash
#SBATCH -J 082025_split_chroms
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 5
#SBATCH -t 06:00:00
#SBATCH --mem=16G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p hotel
#SBATCH -q hotel
#SBATCH -A htl195
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue
# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_SAMTOOLS"

# Get the chromosome name from the command line argument
CHROM=$1

# Your BAM files
BAM1="403_ear_deep_file_hifi_041425_1.bam"
BAM2="403_ear_deep_file_hifi_041425_2.bam"

echo "Starting processing for chromosome: $CHROM"

# Extract reads from BAM1 for the target chromosome
echo "Splitting $BAM1 for $CHROM"
samtools view -b -@ 4 "$BAM1" "$CHROM" > "${BAM1%.bam}.${CHROM}.bam"
samtools index "${BAM1%.bam}.${CHROM}.bam"

# Extract reads from BAM2 for the target chromosome
echo "Splitting $BAM2 for $CHROM"
samtools view -b -@ 4 "$BAM2" "$CHROM" > "${BAM2%.bam}.${CHROM}.bam"
samtools index "${BAM2%.bam}.${CHROM}.bam"

echo "Finished processing for chromosome: $CHROM"
