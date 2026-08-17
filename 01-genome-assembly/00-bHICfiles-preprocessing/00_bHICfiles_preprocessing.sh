#!/bin/bash
#SBATCH -J bhic_preprocess
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=178G
#SBATCH -t 4-00:00:00
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p hotel
#SBATCH -q hotel
#SBATCH -A htl195
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_BULK_HIC"

# Usage: submitted once per sample by
# 00_submittingFiles_for_bHICfiles_preprocessing.sh with three arguments:
#   $1 = read_file1 (raw Hi-C R1), $2 = read_file2 (raw Hi-C R2), $3 = sample name
read_file1="$1"
read_file2="$2"
name1="$3"

genome="octDeg2"                # label only
genome_file="$BHIC_GENOME"
chrom_file="$BHIC_GENOME_FAI"

input_dir="$BHIC_RAW_INPUT_DIR"
data_dir="$BHIC_DATA_DIR"
map_dir_sam="$BHIC_MAP_DIR/$name1"
mat_dir_sam="$BHIC_MATRIX_DIR/$name1"
mkdir -p "$map_dir_sam" "$mat_dir_sam"

read1="${read_file1%.fastq.gz}"
read2="${read_file2%.fastq.gz}"

# (1) Trim adapters/low-quality bases (skip if already done)
if [[ ! -f "$data_dir/${read1}_val_1.fq.gz" && ! -f "$data_dir/${read2}_val_2.fq.gz" ]]; then
  echo "Output file does not exist. Running Trim Galore!"
  trim_galore --cores 8 --paired "$input_dir/$read_file1" "$input_dir/$read_file2" -o "$data_dir"
else
  echo "Output file already exists. Skipping Trim Galore!"
fi

# (2) Index the genome with bwa-mem2 (if not already indexed)
index_files=("${genome_file}.bwt.2bit.64" "${genome_file}.bwt.8bit.32" "${genome_file}.pac")
for file in "${index_files[@]}"; do
  [[ ! -f "$file" ]] && { echo "Building index for bwa-mem2..."; bwa-mem2 index "$genome_file"; break; }
done

# (3) Align and sort
bwa-mem2 mem -SP5M -T0 -t8 "$genome_file" "$data_dir/${read1}_val_1.fq.gz" "$data_dir/${read2}_val_2.fq.gz" \
  | samtools view -bhS - > "$map_dir_sam/${name1}_${genome}.bam"

# (4) Parse pairs, sort, dedup, split
echo "starting pairsam generation"
samtools view -h "$map_dir_sam/${name1}_${genome}.bam" \
  | pairtools parse --min-mapq 40 --walks-policy all --nproc-in 4 --nproc-out 4 \
      --max-inter-align-gap 30 --chroms-path "$chrom_file" --assembly "$genome_file" \
      --output-stats "$map_dir_sam/${name1}_${genome}.pairparse.txt" \
  | pairtools sort --nproc 8 --tmpdir "$SCRATCH_DIR/temp-dir-fast/" \
  > "$map_dir_sam/${name1}_${genome}.sorted.pairsam"

pairtools dedup --nproc-in 4 --nproc-out 4 --mark-dups \
    --output-stats "$map_dir_sam/${name1}_${genome}.pairdedup.txt" \
    "$map_dir_sam/${name1}_${genome}.sorted.pairsam" \
  | pairtools split --nproc-in 4 --nproc-out 4 \
      --output-pairs "$map_dir_sam/${name1}_${genome}.pairs" --output-sam - \
  | samtools view -bS -@ 8 \
  | samtools sort -T "$map_dir_sam" -@ 8 -o "$map_dir_sam/${name1}_${genome}.pairtools.bam"

# Remove the huge intermediate pairsam file
rm "$map_dir_sam/${name1}_${genome}.sorted.pairsam"

# (5) Build and balance contact matrices
bgzip -f "$map_dir_sam/${name1}_${genome}.pairs"
pairix -f "$map_dir_sam/${name1}_${genome}.pairs.gz"

bs=5000
cooler cload pairix "${chrom_file}:${bs}" "$map_dir_sam/${name1}_${genome}.pairs.gz" "$mat_dir_sam/${name1}_${genome}_${bs}.cool"

cooler zoomify --balance --balance-args '--convergence-policy store_nan' -p 8 \
  -o "$mat_dir_sam/${name1}_${genome}.mcool" \
  -r 5000,10000,25000,50000,100000,250000,500000,1000000,2500000 \
  "$mat_dir_sam/${name1}_${genome}_${bs}.cool"
