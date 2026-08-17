#!/bin/bash
#SBATCH -J 100624_HiCAT #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 8 #Number of threads per process
#SBATCH -t 1-00:00:00 #Short for --time walltime limit
#SBATCH --mem=300G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #platinum #Partition name
#SBATCH -q condo #hcp-csd788 #QOS name
#SBATCH -A csd788  #htl195 for hotel #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user you@example.com #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc
conda activate toolshed-HiCAT

hicat --output_dir 300Gb -i /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.chr4.125-135mb.fasta -t /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-centraAnno/hifiasm-0414/chr4_monomerTemplates_line1743-1804.fa -th 8
