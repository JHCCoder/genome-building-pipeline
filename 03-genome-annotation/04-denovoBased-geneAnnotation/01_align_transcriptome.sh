#!/bin/bash
#SBATCH -J 042024_hisatAlign_samtoolSort
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 08:00:00
#SBATCH --mem=64G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user jhc103@ucsd.edu

source /tscc/nfs/home/jhc103/.bashrc
conda activate genome-annotation

# Get SRR ID from job array index
SRR_LIST=("$@")
srr=${SRR_LIST[$SLURM_ARRAY_TASK_ID]}

#may not need if export works: assembly="hifiasm_121624_haphic"
#may not need if export works: base_name="hifiasm_121624_haphic_masked"
assembly_file="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.masked.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked/assembly_final.sorted.headerRenamed.fasta.masked" # Previous assembly file: "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-repeatmasker/${assembly}/${base_name}"
input_mRNA_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/sequencing-mRNAseq/SRA-ncbi-deposits"
output_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-hisat2Aligned-samtoolSorted-mRNA-transcript"
output_dir=${output_dir}/${base_name}
mkdir -p ${output_dir}
mkdir -p ./hisat2_samtools_log
output_log_dir="./hisat2_samtools_log/${base_name}_${srr}"
mkdir -p $output_log_dir

# Run hisat2 alignment
hisat2 -x ${base_name} -1 ${input_mRNA_dir}/${srr}_1.fastq.gz -2 ${input_mRNA_dir}/${srr}_2.fastq.gz --dta -p 12 -S $output_dir/${srr}.sam 1> ${output_log_dir}/hisat2.${base_name}.${srr}.stdout 2> ${output_log_dir}/hisat2.${base_name}.${srr}.stderr

# Run samtools sort
samtools sort -o $output_dir/${srr}.bam -@ 12 $output_dir/${srr}.sam 1> ${output_log_dir}/samtools.${base_name}.${srr}.stdout 2> ${output_log_dir}/samtools.${base_name}.${srr}.stderr
