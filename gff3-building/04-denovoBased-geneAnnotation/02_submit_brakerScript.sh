#!/bin/bash
#SBATCH -J 042425_braker3_full_run #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 48 #Number of threads per process
#SBATCH -t 3-00:00:00 #Short for --time walltime limit
#SBATCH --mem=256G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #Partition name
#SBATCH -q condo #QOS name
#SBATCH -A csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address

module load singularitypro/3.11

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/genome-annotation/braker3
input_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-hisat2Aligned-samtoolSorted-mRNA-transcript/hifiasm_041425_haphic_masked_curated/"  #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-hisat2Aligned-samtoolSorted-mRNA-transcript/hifiasm_022425_haphic_masked/"
SRR_ids="SRR17216301,SRR17216302,SRR17216303,SRR17216304,SRR17216305,SRR17216306,SRR17216307,SRR17216308,SRR17216309,SRR17216310,SRR17216311,SRR17216312,SRR17216313,SRR17216314,SRR17216315,SRR17216293,SRR17216294,SRR17216295,SRR17216296,SRR17216297,SRR17216320,SRR17216319,SRR17216318,SRR17216317,SRR17216316,SRR17216299,SRR17216298,SRR17216321,SRR17216300"
input_SRR_bam=$(echo "$SRR_ids" | sed "s|,|.bam,$input_dir/|g; s|^|$input_dir/|; s|$|.bam|")
wd="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/genome-annotation/braker3/braker-hifiasm-041425-vertebrate-plus-relatives-curated-genome/"
mkdir -p $wd

## Test1 tried with inputting the fastq directly but unparallelized alignments takes way too long!
#singularity exec -B ${PWD}:${PWD} braker3.sif braker.pl --genome=scaffolds.group1.fa.masked --prot_seq=Vertebrata.fa --rnaseq_sets_ids=${SRR_ids} --rnaseq_sets_dirs=${input_dir}/sequencing-mRNAseq/SRA-ncbi-deposits --gff3
## Test2 on one chromosome with aligned bam files
#singularity exec -B ${PWD}:${PWD} braker3.sif braker.pl --genome=scaffolds.group1.fa.masked --prot_seq=Vertebrata.fa --rnaseq_sets_ids=${SRR_ids} --rnaseq_sets_dirs=${input_dir}/sequencing-mRNAseq/SRA-ncbi-deposits --gff3 --threads=48

## Run annotation on vertebrates and relatives!  
singularity exec -B /tscc/lustre/ddn/scratch/jhc103/temp-dir-fast,/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj braker3.sif braker.pl --genome=assembly_final.sorted.headerRenamed.chrAssigned.masked.fasta --prot_seq=Vertebrata_plus_relatives.fa --bam=${input_SRR_bam} --gff3 --threads=48 --workingdir=${wd}

