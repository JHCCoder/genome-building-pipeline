#!/bin/bash
#SBATCH -J 091625_haphic_octDeg1 #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 8 #Number of threads per process
#SBATCH -t 7-00:00:00 #Short for --time walltime limit ==> took less than 3 hours!
#SBATCH --mem=180G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p platinum #Partition name
#SBATCH -q hcp-csd788 #QOS name
#SBATCH -A csd788 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address
#SBATCH --no-requeue

source /tscc/nfs/home/jhc103/.bashrc

ref_file="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/OctDegus1_genome/OctDeg1/fasta/genome.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-verkko/hifi-hic-4thAttempt-moreMemCleanedScratch-48cpus/asm/assembly.haplotype1.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/0414250-assembly/hifiasm-041425-assembly-mitoFiltered/genome_chrom.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/022425-assembly/hifiasm-assembled/deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_022425.asm.hic.p_ctg.fa" ## Genome in december:"/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/data/denovo_OctDegus_genome/121624-assembly/deNovo_hic_ONT_aggressivePurge3_kmer21_121624.hic.p_ctg.fa"
#hic_dir="/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/data/sequencing-reads-HiC"
output_dir="/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/output/outputs-from-haphic-alignment"
hic_file1="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/haphic/input-hic-read/WB_438_R1_trimmed_combined.fq.gz" #"WB_438_2_S1_L008_R1_001_val_1.fq.gz"
hic_file2="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/haphic/input-hic-read/WB_438_R2_trimmed_combined.fq.gz" #"WB_438_2_S1_L008_R2_001_val_2.fq.gz"
prefix="octDeg1" #"references_verkko_male403_hifiHiCMode_081425_hap1"
output_dir=$output_dir/$prefix
mkdir -p $output_dir
cd $output_dir

conda activate bulk-HiC-processing
# (1) Align Hi-C data to the assembly if not exist yet, remove PCR duplicates and filter out secondary and supplementary alignments
if [ -f "${ref_file}.amb" ] && [ -f "${ref_file}.ann" ] && [ -f "${ref_file}.bwt" ] && [ -f "${ref_file}.pac" ] && [ -f "${ref_file}.sa" ]; then
    echo "BWA index files already exist."
else
    echo "BWA index files do not exist. Creating index..."
    bwa index "$ref_file"
fi
#bwa mem -5SP -t 48 $ref_file $hic_dir/$hic_file1 $hic_dir/$hic_file2 | samblaster | samtools view - -@ 48 -S -h -b -F 3340 -o ${output_dir}/${prefix}.bam
bwa mem -5SP -t 8 $ref_file $hic_file1 $hic_file2 | samblaster | samtools view - -@ 8 -S -h -b -F 3340 -o ${output_dir}/${prefix}.bam

conda activate haphic
# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
/tscc/projects/ps-renlab2/jhc103/toolshed/HapHiC/utils/filter_bam ${output_dir}/${prefix}.bam 1 --nm 3 --threads 8 | samtools view - -b -@ 8 -o ${output_dir}/${prefix}.filtered.bam

echo "Finished alignment"

/tscc/projects/ps-renlab2/jhc103/toolshed/HapHiC/haphic pipeline $ref_file ${output_dir}/${prefix}.filtered.bam 30 --threads 4 --processes 2 ## MboI as default--RE "GATC" ## 28 autosomes and 2 sex chromosomes ## Input should be haploid number 
