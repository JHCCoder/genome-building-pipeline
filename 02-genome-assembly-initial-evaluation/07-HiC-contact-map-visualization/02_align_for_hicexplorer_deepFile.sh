#!/bin/bash
#SBATCH -J 080825_align_hicRead_deepDile
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 9
#SBATCH --mem=180G
#SBATCH -t 7-00:00:00
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p hotel
#SBATCH -q hotel
#SBATCH -A htl195
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_BULK_HIC"

genome="hifi_041425"
genome_file="$HIC_GENOME"
chrom_file="${genome_file}.fai"

read1="WB_438_2_S1_L008_R1_001.fq.gz" #"WB_438_R1_trimmed_combined.fq.gz" #"WB_438_2_S6_L003_R1_001.fastq.gz" #"WB_438_2_S1_L008_R1_001_val_1_subsampled.fq.gz" #"WB_438_R1_trimmed_combined.fq.gz"
read2="WB_438_2_S1_L008_R2_001.fq.gz" #"WB_438_R2_trimmed_combined.fq.gz" #"WB_438_2_S6_L003_R2_001.fastq.gz" #"WB_438_2_S1_L008_R2_001_val_2_subsampled.fq.gz" #"WB_438_R2_trimmed_combined.fq.gz"
run="403_ear_deep_file"
data_dir="$DATA_DIR/sequencing-reads-HiC" #"$HIC_HIFI_READ_DIR" #"$DATA_DIR/sequencing-reads-HiC" #"$HIC_HIFI_READ_DIR"
output_dir="$HIC_OUT_DIR"

echo "read_file1: $read1"
echo "read_file2: $read2"
name1=$run

# Generate base names without .fq.gz extension
read1_base="${read1%.fq.gz}"
read2_base="${read2%.fq.gz}"

# Check for trimmed files
if [ ! -f "${data_dir}/${read1_base}_val_1.fq.gz" ] || [ ! -f "${data_dir}/${read2_base}_val_2.fq.gz" ]; then
    echo "Running Trim Galore..."
    trim_galore --cores 8 --paired "${data_dir}/${read1}" "${data_dir}/${read2}" -o "${data_dir}"
else
    echo "Trimmed files exist. Skipping Trim Galore!"
fi

# Check for BWA-MEM2 index
if [ ! -f "${genome_file}.bwt.2bit.64" ]; then
    echo "Building index for bwa-mem2..."
    bwa-mem2 index "${genome_file}"
fi

# Ensure output directory exists
mkdir -p "${output_dir}"

# Align reads
echo "Aligning read 1..."
bwa-mem2 mem -A 1 -B 4 -E 50 -L 0 -t 8 "${genome_file}" "${data_dir}/${read1_base}_val_1.fq.gz" | samtools view -Shb - | \
samtools sort -@ 2 -m 2G -T "${output_dir}/tmp_1" -o "${output_dir}/${name1}_${genome}_1.bam"

echo "Aligning read 2..."
bwa-mem2 mem -A 1 -B 4 -E 50 -L 0 -t 8 "${genome_file}" "${data_dir}/${read2_base}_val_2.fq.gz" | samtools view -Shb - | \
samtools sort -@ 2 -m 2G -T "${output_dir}/tmp_2" -o "${output_dir}/${name1}_${genome}_2.bam"

echo "Alignment completed."

# Index the sorted BAM files
echo "Indexing sorted BAM files..."
samtools index -@ 2 "${output_dir}/${name1}_${genome}_1.sorted.bam"
samtools index -@ 2 "${output_dir}/${name1}_${genome}_2.sorted.bam"
