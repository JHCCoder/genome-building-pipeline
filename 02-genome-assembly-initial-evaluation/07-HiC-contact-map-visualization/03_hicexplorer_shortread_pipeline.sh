#!/bin/bash
#SBATCH -J 052125_combined_hic_pipeline
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 9
#SBATCH --mem=192G
#SBATCH -t 2-00:00:00
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p platinum
#SBATCH -q hcp-csd788
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load environment and activate conda
# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)

# Configuration
genome="hifi_041425"
genome_file="$HIC_GENOME"
restriction_cut_file="$HIC_RESTRICTION_BED"
runs=("SRR19145623" "SRR19145629")
output_dir="$HIC_OUT_DIR"
input_dir="$HIC_SHORTREAD_DIR"
hic_matrix_dir="hicMatrix"

# Create output directories
mkdir -p "${output_dir}" "${hic_matrix_dir}"

# Step 1: Alignment (using bulk-HiC-processing environment)
conda activate "$ENV_BULK_HIC"

# Check for BWA-MEM2 index
if [ ! -f "${genome_file}.bwt.2bit.64" ]; then
    echo "Building index for bwa-mem2..."
    bwa-mem2 index "${genome_file}"
fi

for run in "${runs[@]}"; do
    # Define paths
    data_dir="${input_dir}/${run}_24threads"
    read1="${run}_1.fastq.gz"
    read2="${run}_2.fastq.gz"
    read1_base="${run}_1"
    read2_base="${run}_2"

    # Trim Galore (if not already done)
    if [ ! -f "${data_dir}/${read1_base}_val_1.fq.gz" ] || [ ! -f "${data_dir}/${read2_base}_val_2.fq.gz" ]; then
        echo "Running Trim Galore for ${run}..."
        trim_galore --cores 8 --paired "${data_dir}/${read1}" "${data_dir}/${read2}" -o "${data_dir}"
    else
        echo "Trimmed files exist for ${run}. Skipping Trim Galore!"
    fi

    # Align reads with read groups
    echo "Aligning ${run} read pairs..."
    bwa-mem2 mem -A 1 -B 4 -E 50 -L 0 -t 8 -R "@RG\tID:${run}_1\tSM:${run}" \
        "${genome_file}" \
        "${data_dir}/${read1_base}_val_1.fq.gz" \
        | samtools view -Shb - > "${output_dir}/${run}_${genome}_1.bam"

    bwa-mem2 mem -A 1 -B 4 -E 50 -L 0 -t 8 -R "@RG\tID:${run}_2\tSM:${run}" \
        "${genome_file}" \
        "${data_dir}/${read2_base}_val_2.fq.gz" \
        | samtools view -Shb - > "${output_dir}/${run}_${genome}_2.bam"

    # Sort and index
    for pair in 1 2; do
        samtools sort -@ 8 -o "${output_dir}/${run}_${genome}_${pair}_sorted.bam" \
            "${output_dir}/${run}_${genome}_${pair}.bam"
        samtools index "${output_dir}/${run}_${genome}_${pair}_sorted.bam"
    done
done

# Step 2: Hi-C matrix building (using toolshed-HiCExplorer environment)
conda activate "$ENV_HICEXPLORER"

for run in "${runs[@]}"; do
    echo "Building Hi-C matrix for ${run}..."
    hicBuildMatrix \
        --samFiles "${output_dir}/${run}_${genome}_1_sorted.bam" \
                   "${output_dir}/${run}_${genome}_2_sorted.bam" \
        --binSize 10000 \
        --restrictionSequence GATC \
        --danglingSequence GATC \
        --restrictionCutFile "${restriction_cut_file}" \
        --outBam "${genome}_ref.bam" \
        --outFileName "${hic_matrix_dir}/${run}_${genome}_10kb.h5" \
        --QCfolder "${hic_matrix_dir}/${run}_${genome}_10kb_QC" \
        --threads 8 \
        --inputBufferSize 400000
done

# Step 3: Merge matrices (optional)
echo "Merging Hi-C matrices..."
hicMerge \
    --matrix "${hic_matrix_dir}/SRR19145623_${genome}_10kb.h5" \
    --matrix "${hic_matrix_dir}/SRR19145629_${genome}_10kb.h5" \
    --outFileName "${hic_matrix_dir}/merged_${genome}_10kb.h5"

echo "Pipeline completed successfully!"
