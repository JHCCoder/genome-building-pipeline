#!/bin/bash
#SBATCH -J 090325_run_trf_single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2-00:00:00       # Increased time since it's processing the whole file
#SBATCH --mem=64G         # Increased memory since it's processing the whole file
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p hotel #platinum #condo
#SBATCH -q hotel #hcp-csd788 #condo
#SBATCH -A htl195 #csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

# Load required modules
module load shared
module load gpu/0.17.3
module load parallel/20210922-hzmc3me
module load pyfasta/0.5.2

# Set variables
INFILE="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/output/outputs-from-haphic-alignment/references_verkko_male403_hifiHiCMode_081425_hap1/04.build/scaffolds.fa" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/genome-related-species/GCA_964261345.1_mHetGlaV3_genomic.fna" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/genome-related-species/GCF_034190915.1_mCavPor4.1_genomic.fna" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/genome-related-species/GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.fna"  #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/genome-related-species/GCF_036323735.1_GRCr8_genomic.fna" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/genome-related-species/GCA_944319725.1_Naked_mole-rat_paternal_genomic.fna" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/GRCm39_genome/GCF_000001635.27_GRCm39_genomic.fna" #"/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/Homo_sapien_genome/hg38/refseq-convention/GCF_000001405.26_GRCh38_genomic.fna"  # Now includes .fasta extension
OUTPUT_DIR="/tscc/lustre/ddn/scratch/jhc103/trf-output-deguVerkko"
mkdir -p $OUTPUT_DIR

# Run TRF on the complete file
echo "Running TRF on complete file..."
cd $OUTPUT_DIR
/tscc/projects/ps-renlab2/jhc103/toolshed/TRF-4.09.1/build/src/trf \
    $INFILE \
    2 7 7 80 10 50 500 -d -l 20

echo "TRF processing complete. Final outputs:"
echo "- Masked sequences: $OUTPUT_DIR/${INFILE}.mask"
echo "- Repeat annotations: $OUTPUT_DIR/${INFILE}.dat"
