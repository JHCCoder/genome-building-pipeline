#!/bin/bash
#SBATCH -J 100225_repeatMasker_verkko #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 33 #Number of threads per process
#SBATCH -t 3-00:00:00 #Short for --time walltime limit
#SBATCH --mem 196G # Request 80 GB of memory
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #Partition name
#SBATCH -q condo #QOS name
#SBATCH -A csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address

source /tscc/nfs/home/jhc103/.bashrc
conda activate toolshed-repeatmodeler

## Make sure the header is less than 50 chars long
input_fa="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_verkko_male403_hifiHiCMode_081425_hap1/04.build/scaffolds.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/octDeg1/04.build/scaffolds.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-purged/purged.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-verkko/hifi-hic-4thAttempt-moreMemCleanedScratch-48cpus/asm/assembly.haplotype1.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded2/assembly_final.sorted.headerRenamed.fasta" #'/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded1/scaffolds.fa' #'/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_hifiasm_male403_hifiHiCMode_022425/04.build/scaffolds.fa'
filename=$(basename "$input_fa")
genome_fa_name="verkko_hap1_scaffolded"
db_name="odegus"

output_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-repeatmasker/${genome_fa_name}"
mkdir -p $output_dir
cd $output_dir
cp $input_fa .

## Step 1: Build database (couple of minutes)
BuildDatabase -name $db_name ./$filename

# Step 2: Build de novo repeat library using the masked assembly (17-18 hours)
RepeatModeler -database $db_name -threads 32 -LTRStruct -engine ncbi > out.log
## Memory Efficiency: 8.80% of 660.00 GB ==> could lower this by alot! 

ln -s RM_*/consensi.fa.classified ./

# Step 3: RepeatMasker: Use de novo repeat library to mask our genome 
RepeatMasker -e rmblast -pa 8 -gff -lib consensi.fa.classified -xsmall ./$filename 
