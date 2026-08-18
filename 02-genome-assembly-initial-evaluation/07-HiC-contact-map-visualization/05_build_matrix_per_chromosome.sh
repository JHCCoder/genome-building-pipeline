#!/bin/bash
#SBATCH -J 082625_hicBuildMatrix
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=178G
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

genome="hifi_041425"
run="403_ear_deep_file"
TMP_ROOT="${TMPDIR:-$WORKDIR/tmp}"
SAMBAMBA_TMP_DIR="$TMP_ROOT/sambamba_${SLURM_JOB_ID}"
# Define chromosomes to process (adjust for your genome)
CHROMOSOMES=("chr1" "chr2")

mkdir -p hicMatrix compartments
mem_gb=$(( SLURM_MEM_PER_NODE / 1000 )) - 2
for chrom in "${CHROMOSOMES[@]}"; do
    echo "Processing $chrom..."
    conda activate "$ENV_SAMTOOLS"
    # 1. Subset BAM files for this chromosome
    #samtools view -b ${run}_${genome}_1.sorted.bam "$chrom" > ${run}_${genome}_1.sorted.${chrom}.bam
    samtools view -b ${run}_${genome}_1.sorted.bam "$chrom" | \
	    sambamba sort -t $SLURM_CPUS_PER_TASK -m ${mem_gb}GB --tmpdir="$SAMBAMBA_TMP_DIR" -N -o ${run}_${genome}_1.${chrom}.nameSorted.bam /dev/stdin

    #samtools view -b ${run}_${genome}_2.sorted.bam "$chrom" > ${run}_${genome}_2.sorted.${chrom}.bam
    samtools view -b ${run}_${genome}_2.sorted.bam "$chrom" | \
	    sambamba sort -t $SLURM_CPUS_PER_TASK -m ${mem_gb}GB --tmpdir="$SAMBAMBA_TMP_DIR" -N -o ${run}_${genome}_2.${chrom}.nameSorted.bam /dev/stdin
    #samtools index ${run}_${genome}_2.${chrom}.nameSorted.bam 
    #samtools index ${run}_${genome}_2.${chrom}.nameSorted.bam

    conda activate "$ENV_HICEXPLORER"
    # 2. Build matrix for this chromosome only
    hicBuildMatrix --samFiles ${run}_${genome}_1.${chrom}.nameSorted.bam ${run}_${genome}_2.${chrom}.nameSorted.bam \
                   --binSize 100000 \
                   --restrictionSequence GATC \
                   --danglingSequence GATC \
                   --restrictionCutFile $HIC_RESTRICTION_BED \
                   --outFileName hicMatrix/${run}_${genome}_100kb_${chrom}.h5 \
                   --QCfolder hicMatrix/QC_${chrom} \
                   --threads 2 \
                   --inputBufferSize 200000

    # 3. Correct matrix
    hicCorrectMatrix correct --matrix hicMatrix/${run}_${genome}_100kb_${chrom}.h5 \
                             --correctionMethod ICE \
                             --outFileName hicMatrix/${run}_${genome}_100kb_${chrom}_corrected.h5

    # 4. Call A/B compartments
    hicPCA -m hicMatrix/${run}_${genome}_100kb_${chrom}_corrected.h5 \
           -o compartments/${run}_${genome}_compartments_${chrom}.bed \
           --format bed \
           --which eigensystem
done

# 5. Combine all chromosome results
cat compartments/${run}_${genome}_compartments_*.bed > compartments/${run}_${genome}_compartments_genomewide.bed
