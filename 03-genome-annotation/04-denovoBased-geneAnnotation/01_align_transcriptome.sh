#!/bin/bash
#SBATCH -J hisat2_align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 08:00:00
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

# Array job: SRR IDs are passed as positional args by the submit script.
SRR_LIST=("$@")
srr="${SRR_LIST[$SLURM_ARRAY_TASK_ID]}"

output_dir="$HISAT2_OUT_DIR/$HISAT2_INDEX_NAME"
log_dir="./hisat2_samtools_log/${HISAT2_INDEX_NAME}_${srr}"
mkdir -p "$output_dir" "$log_dir"

# Align RNA-seq reads to the masked assembly (HISAT2 index built by the submit
# script) and coordinate-sort the output.
hisat2 -x "$HISAT2_INDEX_NAME" \
  -1 "$MRNA_DIR/${srr}_1.fastq.gz" -2 "$MRNA_DIR/${srr}_2.fastq.gz" \
  --dta -p 12 -S "$output_dir/${srr}.sam" \
  1> "${log_dir}/hisat2.${HISAT2_INDEX_NAME}.${srr}.stdout" \
  2> "${log_dir}/hisat2.${HISAT2_INDEX_NAME}.${srr}.stderr"

samtools sort -o "$output_dir/${srr}.bam" -@ 12 "$output_dir/${srr}.sam" \
  1> "${log_dir}/samtools.${HISAT2_INDEX_NAME}.${srr}.stdout" \
  2> "${log_dir}/samtools.${HISAT2_INDEX_NAME}.${srr}.stderr"
