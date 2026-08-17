#!/bin/bash
#SBATCH -J 050825_blastp_run_aa_vertebrate_UNIPROT_1thread_2G_platinum #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 1 #Number of threads per process
#SBATCH -t 1:00:00 #Short for --time walltime limit
#SBATCH --mem=2G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p platinum #condo #Partition name
#SBATCH -q hcp-csd788 #condo #QOS name
#SBATCH -A csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc
#conda activate genome-annotation

#space=`df -T /scratch | grep scratch | awk '{print $5}'`  # I won't add "-h" option to df to make it human-readable otherwise the unit will either show "G" or "T" based on the actual available space
#threshold=$((888 * 1024 * 1024))  # Convert 888 GB to 1K blocks (888 * 1024 * 1024)

: '
if [ "$space" -gt "$threshold" ]; then
        TIC=`date +%s`
        cp -r /cm/shared/apps/data/blast/FASTA/nr /scratch/$USER/job_${SLURM_JOBID}  
        cp -r /cm/shared/apps/data/blast/FASTA/taxdb.bti /scratch/$USER/job_${SLURM_JOBID} 
        cp -r /cm/shared/apps/data/blast/FASTA/taxdb.btd /scratch/$USER/job_${SLURM_JOBID} 
        cp -r /cm/shared/apps/data/blast/FASTA/taxonomy4blast.sqlite3 /scratch/$USER/job_${SLURM_JOBID} 
        TOC=`date +%s`
        J1_TIME=$((TOC - TIC))
        echo "transfer took=$J1_TIME sec"
fi


## Copy to project directory if not there 
if [ ! -d "/tscc/projects/ps-renlab2/share/ncbi-nr-database/nr" ]; then
	TIC=`date +%s`
	cp -r /cm/shared/apps/data/blast/FASTA/nr /tscc/projects/ps-renlab2/share/ncbi-nr-database 
	cp -r /cm/shared/apps/data/blast/FASTA/taxdb.bti /tscc/projects/ps-renlab2/share/ncbi-nr-database
	cp -r /cm/shared/apps/data/blast/FASTA/taxdb.btd /tscc/projects/ps-renlab2/share/ncbi-nr-database
	cp -r /cm/shared/apps/data/blast/FASTA/taxonomy4blast.sqlite3 /tscc/projects/ps-renlab2/share/ncbi-nr-database
	TOC=`date +%s`
	J1_TIME=$((TOC - TIC))
	echo "transfer took=$J1_TIME sec"
fi
'
echo "running blastp"
export BLASTDB=/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/data/protein-databases/uniprot-vertebrate #/tscc/projects/ps-renlab2/share/ncbi-nr-database/nr_mammalian_db #/tscc/projects/ps-renlab2/share/ncbi-nr-database/nr #/cm/shared/apps/data/blast/FASTA/nr #/tscc/projects/ps-renlab2/share/ncbi-nr-database/nr  # lustre scratch copy: /tscc/lustre/ddn/scratch/jhc103/temp-dir-fast/nr-database/nr #yuwu's copy: /tscc/lustre/ddn/scratch/ychen64/blast/FASTA/nr 
/tscc/projects/ps-renlab2/jhc103/toolshed/blast+/ncbi-blast-2.16.0+/bin/blastp -query braker_121624_100_entries.aa -db UniProtVertebrate_DB -out /tscc/lustre/ddn/scratch/jhc103/output-blastp-test/braker_121624_100_entries_vertebrateUniprot_database_thread1.tsv -evalue 1e-5 -num_threads 1 -max_target_seqs 8 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -seg yes
#blastp -query ../braker3/braker-121624-output/braker.aa -db nr -out blastp_results.tsv -evalue 1e-5 -num_threads 16 -max_target_seqs 10 -outfmt "6 qaccver saccver scomname pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -seg yes
