#!/bin/bash
#SBATCH -J 071425_sry_blast #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 1 #Number of threads per process
#SBATCH -t 01:00:00 #Short for --time walltime limit
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
assembly_alias="hifiasm_041425_scaffolded_curated_masked"
assembly='/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked/assembly_final.sorted.headerRenamed.fasta.masked' #'/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded2/assembly_final.sorted.fasta' #'/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly/deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_041325.asm.hic.p_ctg.fa' #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_hifiasm_male403_hifiHiCMode_041425_trail1_allChrom/04.build/scaffolds.fa"  #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/022425-assembly/hifiasm-assembled/deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_022425.asm.hic.p_ctg.fa"
sry_gene="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/assigning-chromosomes/sry-gene/all_sry_gene.fasta" #'/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/assigning-chromosomes/sry-gene/rodent_sry_gene.fasta' #'/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/assigning-chromosomes/sry-gene/XM_010621263.1_damara_mole_rat.fasta' #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/assigning-chromosomes/y-chromosomes/Ychrom.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/degu_mitochondria/octDeg1_mitochondria_chromosome.fasta"
# If we are using another assembly! 
## ==> makeblastdb -in $assembly -dbtype nucl
#blastn -query $sry_gene -db $assembly -outfmt "6 qseqid sseqid pident qcovs qcovhsp length gapopen evalue bitscore mismatch qstart qend sstart send" -evalue 1e-5 > ${assembly_alias}_all_sry_gene.out ##  The default output format is -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
blastn -query $sry_gene -db $assembly -outfmt 5 -out results_for_kablammo.xml



