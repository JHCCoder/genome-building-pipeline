#!/bin/bash
#SBATCH -J 090725_multiSubmission_bulk_HiC_preprocessing_403_ear_deep_2_octDeg1_continued #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 8 #Number of threads per process
#SBATCH --mem=178G
#SBATCH -t 4-00:00:00 #Short for --time walltime limit
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p hotel #Partition name
#SBATCH -q hotel #QOS name
#SBATCH -A htl195 #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user jhc103@ucsd.edu #Optional, Send mail to this address

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

genome="octDeg2" #Other genome tested: "hifiasm_121624"
genome_file="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/OctDegus1_genome/OctDeg1/fasta/genome.fa" #Path to other genome: "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/assembly-of-interest/hifiasm-121624/deNovo_hic_ONT_aggressivePurge3_kmer21_121624.hic.p_ctg.fa"
chrom_file="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta.fai" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/OctDegus1_genome/OctDeg1/fasta/genome.fa.fai" # Path to other genome chrom file: "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/assembly-of-interest/hifiasm-121624/deNovo_hic_ONT_aggressivePurge3_kmer21_121624.hic.p_ctg.fa.fai"


map_dir="/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/output/mapped_alignments"

mat_dir="/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/output/matrix"
input_dir="/tscc/projects/ps-renlab2/nzemke/Element_share/20250811"
data_dir="/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/data/sequencing-reads-HiC"

echo "read_file1: $1"
echo "read_file2: $2"
## This is the name of the sample 
name1=$3 ##$(echo "$read_file1" | awk -F '_' '{print $1 "_" $2}')
## These are name of reads 
read1=${1%.fastq.gz}
read2=${2%.fastq.gz}

map_dir_sam=$map_dir/$name1
mat_dir_sam=$mat_dir/$name1
mkdir -p $map_dir_sam
#mkdir -p $trim_dir_sam
mkdir -p $mat_dir_sam

# Check if the output trim file exists
if [ ! -f "$data_dir/${read1}_val_1.fq.gz" ] && [ ! -f "$data_dir/${read1}_val_1.fq.gz" ]; then
    echo "Output file does not exist. Running Trim Galore!"
    trim_galore --cores 8 --paired $input_dir/$1 $input_dir/$2 -o $data_dir
else
    echo "Output file already exists. Skipping Trim Galore!"
fi

## Check if the index files are available
index_files=("${genome_file}.bwt.2bit.64" "${genome_file}.bwt.8bit.32" "${genome_file}.pac")
for file in "${index_files[@]}"; do
    [[ ! -f "$file" ]] && { echo "Building index for bwa-mem2..."; bwa-mem2 index "$genome_file"; break; }
done

## Generate alignment then pair sorted alignment
bwa-mem2 mem -SP5M -T0 -t8 $genome_file $data_dir/${read1}_val_1.fq.gz $data_dir/${read2}_val_2.fq.gz | samtools view -bhS - > $map_dir_sam/${name1}_${genome}.bam


echo "starting pairsam generation"
samtools view -h $map_dir_sam/${name1}_${genome}.bam | pairtools parse --min-mapq 40 --walks-policy all --nproc-in 4 --nproc-out 4 --max-inter-align-gap 30 --chroms-path $chrom_file --assembly $genome_file --output-stats ${map_dir_sam}/${name1}_${genome}.pairparse.txt | pairtools sort --nproc 8 --tmpdir /tscc/lustre/ddn/scratch/jhc103/temp-dir-fast/ > ${map_dir_sam}/${name1}_${genome}.sorted.pairsam

pairtools dedup --nproc-in 4 --nproc-out 4 --mark-dups --output-stats ${map_dir_sam}/${name1}_${genome}.pairdedup.txt ${map_dir_sam}/${name1}_${genome}.sorted.pairsam | pairtools split --nproc-in 4 --nproc-out 4 --output-pairs ${map_dir_sam}/${name1}_${genome}.pairs --output-sam - | samtools view -bS -@ 8 | samtools sort -T ${map_dir_sam} -@ 8 -o ${map_dir_sam}/${name1}_${genome}.pairtools.bam

###### 
## Delete the huge pairsam file 
rm ${map_dir_sam}/${name1}_${genome}.sorted.pairsam

################
## Generate matrix files 

bgzip -f ${map_dir_sam}/${name1}_${genome}.pairs
pairix -f ${map_dir_sam}/${name1}_${genome}.pairs.gz

# ### binning
bs=5000
#cooler cload pairix ${chrom_file}:${bs} ${map_dir_sam}/${name1}_${genome}.pairs.gz ${mat_dir_sam}/${name1}_${genome}_${bs}.cool

# ### Balancing
cooler zoomify --balance --balance-args '--convergence-policy store_nan' -p 8 -o ${mat_dir_sam}/${name1}_${genome}.mcool -r 5000,10000,25000,50000,100000,250000,500000,1000000,2500000 ${mat_dir_sam}/${name1}_${genome}_${bs}.cool

