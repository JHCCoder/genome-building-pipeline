#!/bin/bash
#SBATCH -J 102424_genome_annotation_transfer #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 4 #Number of threads per process
#SBATCH -t 01:00:00 #Short for --time walltime limit
#SBATCH --mem=80G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #Partition name
#SBATCH -q condo #QOS name
#SBATCH -A csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address

source /tscc/nfs/home/jhc103/.bashrc
conda activate toolshed-liftoff

assembly="verkko-hap1-transferred-from-hifiasm-041425" #"hifiasm-041425-scaffolded-chrAssigned-mito"
query_fa="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-repeatmasker/verkko_hap1_scaffolded/scaffolds.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded2/assembly_final.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/0414250-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded1/scaffolds.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/022425-assembly/hifiasm-assembled-scaffolded-masked/scaffolds.fa.masked" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-repeatmasker/$assembly/scaffolds.fa.masked"
ref_fa="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/OctDegus1_genome/OctDeg1/fasta/GCF_000260255.1_OctDeg1.0_genomic.fna"
ref_gff="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm_041425_denovoEnhanced_peaks2utr_sorted.gff3" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/OctDegus1_genome/OctDeg1/genes/GCF_000260255.1_OctDeg1.0_genomic.gff"
output_dir="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-liftoff/$assembly"
mkdir -p $output_dir
temp_dir="/tscc/lustre/ddn/scratch/jhc103/temp-dir-fast/liftoff-intermediate-files/$assembly"
mkdir -p $temp_dir

liftoff $query_fa $ref_fa -g $ref_gff -o $output_dir/$assembly.gff -u $output_dir/unmapped_features.txt -dir $temp_dir -p 4 > liftoff_${assembly}_run.out
