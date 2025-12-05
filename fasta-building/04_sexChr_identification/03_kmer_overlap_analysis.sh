#!/bin/bash
#SBATCH -J 102425_discoverY #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 2 #Number of threads per process
#SBATCH -t 12:00:00 #Short for --time walltime limit
#SBATCH --mem=360G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #Partition name
#SBATCH -q condo #QOS name
#SBATCH -A csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc
conda activate toolshed-DiscoverY

## Used the bloom filter
#python discoverY.py --female_kmers_set --female_bloom --kmer_size 25 --mode female+male --female_bloom_capacity 3100000000

python discoverY.py --female_kmers_set --kmer_size 25 --mode female+male
