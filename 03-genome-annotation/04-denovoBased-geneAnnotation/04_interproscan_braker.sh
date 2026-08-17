#!/bin/bash
#SBATCH -J 050925_interproscan_stragglers #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 4 #Number of threads per process
#SBATCH -t 06:00:00 #Short for --time walltime limit
#SBATCH --mem=30G
#SBATCH --array=66     # 10 array jobs, max 4 running at once
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #platinum #Partition name
#SBATCH -q condo #hcp-csd788 #QOS name
#SBATCH -A csd788 #htl195 for hotel #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc
#conda activate r-seurat


sed 's/\*//g' 041425-braker-split-files-100files/chunk_${SLURM_ARRAY_TASK_ID}.fa > 041425-braker-interproscan-result/chunk_${SLURM_ARRAY_TASK_ID}.clean.aa
/tscc/projects/ps-renlab2/jhc103/toolshed/interproscan/interproscan-5.74-105.0/interproscan.sh -i 041425-braker-interproscan-result/chunk_${SLURM_ARRAY_TASK_ID}.clean.aa -f tsv -dp -o 041425-braker-interproscan-result/interproscan_result_${SLURM_ARRAY_TASK_ID} -cpu 4
