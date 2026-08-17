#!/bin/bash
#SBATCH -J mito_chrom_blast #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 1 #Number of threads per process
#SBATCH -t 03:00:00 #Short for --time walltime limit
#SBATCH --mem=64G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #Partition name
#SBATCH -q condo #QOS name
#SBATCH -A csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address

source /tscc/nfs/home/jhc103/.bashrc
conda activate genome-annotation
assembly_alias="hifiasm_041425"
assembly="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/0414250-assembly/hifiasm-041425-assembly/deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_041325.asm.hic.p_ctg.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/022425-assembly/hifiasm-assembled/deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_022425.asm.hic.p_ctg.fa"
mito_chrom="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/degu_mitochondria/octDeg1_mitochondria_chromosome.fasta"
makeblastdb -in $assembly -dbtype nucl
blastn -query $mito_chrom -db $assembly -outfmt "6 qseqid sseqid pident qcovs length gapopen evalue bitscore mismatch" -evalue 1e-5 > ${assembly_alias}_mito_blast.out ##  The default output format is -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"




