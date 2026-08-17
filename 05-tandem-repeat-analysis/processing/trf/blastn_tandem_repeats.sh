#!/bin/bash
#SBATCH -J 072025_blastn_all_seq #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 4 #Number of threads per process
#SBATCH -t 7-00:00:00 #Short for --time walltime limit
#SBATCH --mem=128G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p hotel #condo #platinum #Partition name
#SBATCH -q hotel #condo #hcp-csd788 #QOS name
#SBATCH -A htl195 #csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user you@example.com #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc
conda activate genome-annotation

# Configuration
input_file="389peak_repeat.fasta" # "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/tandem-repeat-investigation/degu_peaks_length100_periodsize190_repeat.fasta" #"389peak_repeat.fasta"
backup_database="/cm/shared/apps/data/blast/FASTA/nt"
OUTDIR="./blastn-results"
mkdir -p "$OUTDIR"
out_name=$(basename "$input_file" .fasta)


export BLASTDB="$backup_database"
#blastn -query $input_file -db nt -dust no  -out $OUTDIR/$out_name.txt -outfmt 6 -evalue 1e-5
blastn -query $input_file -db nt -dust no -word_size 7 -evalue 1e-5 -outfmt 6 -out $OUTDIR/${out_name}_cross_species_hits_more_sequence.txt -num_threads "$SLURM_CPUS_PER_TASK"
