#!/bin/bash
#SBATCH -J 062425_cd_hit #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 32 #Number of threads per process
#SBATCH -t 24:00:00 #Short for --time walltime limit
#SBATCH --mem=64G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #platinum #Partition name
#SBATCH -q condo #hcp-csd788 #QOS name
#SBATCH -A csd788  #htl195 for hotel #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user you@example.com #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc
conda activate genome-annotation

cd-hit-est -i assembly_final.sorted.headerRenamed.chrAssigned.monomers.fa -o assembly_final.sorted.headerRenamed.chrAssigned.monomers.clustered.fa -c 0.9 -n 5 -T 32 -M 64000

echo finished
