#!/bin/bash
#SBATCH -J 092525_centAnno #Optional, short for --job-name
#SBATCH --nodelist=tscc-14-12
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 8 #Number of threads per process
#SBATCH -t 1-00:00:00 #Short for --time walltime limit
#SBATCH --mem=48G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p platinum #Partition name
#SBATCH -q hcp-csd788 #QOS name
#SBATCH -A csd788  #htl195 for hotel #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user you@example.com #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc

./centroAnno /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/Homo_sapien_genome/hg38/refseq-convention/GCF_000001405.26_GRCh38_genomic.fna -o ~/ps-renlab2-link/degu-genome-assembly-proj/output/outputs-from-centraAnno/hg38 -x anno-asm -t 8
