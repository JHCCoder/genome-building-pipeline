#!/bin/bash
#SBATCH -J 020324_repeatModeler #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 33 #Number of threads per process
#SBATCH -t 7-00:00:00 #Short for --time walltime limit
#SBATCH --mem 660G # Request 80 GB of memory
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #Partition name
#SBATCH -q condo #QOS name
#SBATCH -A csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address

source /tscc/nfs/home/jhc103/.bashrc
conda activate toolshed-repeatmodeler

## Build database
BuildDatabase -name odegus /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_hifiasm_male403_hifiHiCMode_022425/04.build/scaffolds.fa

# Step 2: Build de novo repeat library using the masked assembly
# Use the soft-masked output from step 1
initial_masked_genome="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_hifiasm_male403_hic_ont_121624/04.build/scaffolds.fa"
#temp_dir="/tscc/lustre/ddn/scratch/jhc103/repeatModel-odegus"
output_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-repeatmasker/hifiasm-121624-haphic"

mkdir -p $output_dir
#cd $output_dir

#BuildDatabase -engine ncbi -name odegus_db $initial_masked_genome 
RepeatModeler -database odegus -threads 32 -LTRStruct -engine ncbi > out.log

## Move result out of temp
cp -r odegus_db* $output_dir
