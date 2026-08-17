#!/bin/bash
#SBATCH -J 052725_hifiasmHiCMode_hifiData_lightPurge1_kmer21
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH -t 4-00:00:00
#SBATCH --mem=256G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p platinum #condo
#SBATCH -q hcp-csd788 #condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user jhc103@ucsd.edu

output_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-hifiasm"
file_name="deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_052325.asm" #"deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_022425.asm" #"deNovo_hifiHiCMode_hifiData4files_lightPurge1_kmer21_022425.asm"

output_file="${output_dir}/${file_name}"

/tscc/projects/ps-renlab2/jhc103/toolshed/hifiasm/hifiasm -o $output_file --primary -t36 -l3 -k21 --h1 /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/sequencing-reads-HiC/WB_438_2_S1_L008_R1_001_val_1.fq.gz --h2 /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/sequencing-reads-HiC/WB_438_2_S1_L008_R2_001_val_2.fq.gz /tscc/lustre/ddn/scratch/jhc103/hifi/male/*/*fastq # This is the other file with the ONT sequenced reads: /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-dorado-basecaller/all_male458_samples.fq.gz 


