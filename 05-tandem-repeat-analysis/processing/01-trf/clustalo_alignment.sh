#!/bin/bash
#SBATCH -J 073025_clustalo_aligned #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 8 #Number of threads per process
#SBATCH -t 24:00:00 #Short for --time walltime limit
#SBATCH --mem=128G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p hotel #condo #platinum #Partition name
#SBATCH -q hotel #condo #hcp-csd788 #QOS name
#SBATCH -A htl195 #csd788  #htl195 for hotel #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user you@example.com #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc
conda activate evolutionary-tree

clustalo --threads=8 -t DNA -i 349peak_repeat_10Longest.fasta -o 349peak_repeat_10Longest_aligned.fasta --outfmt=fasta --force
echo testing
