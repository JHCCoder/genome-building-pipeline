#!/bin/bash
#SBATCH -J haphic_scaffold
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 7-00:00:00
#SBATCH --mem=180G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p platinum
#SBATCH -q hcp-csd788
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)

ref_file="$HAPHIC_REF"     # assembly to scaffold
hic_file1="$HIC_HAPHIC_R1"
hic_file2="$HIC_HAPHIC_R2"
prefix="octDeg1"
output_dir="$OUTPUT_DIR/outputs-from-haphic-alignment/$prefix"
mkdir -p "$output_dir"
cd "$output_dir"

conda activate "$ENV_BULK_HIC"

# (1) Align Hi-C reads to the assembly (index first if needed), remove PCR
#     duplicates, and drop secondary/supplementary alignments.
if [[ -f "${ref_file}.amb" && -f "${ref_file}.ann" && -f "${ref_file}.bwt" \
      && -f "${ref_file}.pac" && -f "${ref_file}.sa" ]]; then
  echo "BWA index files already exist."
else
  echo "BWA index files do not exist. Creating index..."
  bwa index "$ref_file"
fi
bwa mem -5SP -t 8 "$ref_file" "$hic_file1" "$hic_file2" | samblaster | \
  samtools view - -@ 8 -S -h -b -F 3340 -o "${output_dir}/${prefix}.bam"

conda activate "$ENV_HAPHIC"

# (2) Filter alignments: keep MAPQ >= 1 and edit distance < 3.
"$HAPHIC_DIR/utils/filter_bam" "${output_dir}/${prefix}.bam" 1 --nm 3 --threads 8 | \
  samtools view - -b -@ 8 -o "${output_dir}/${prefix}.filtered.bam"

echo "Finished alignment"

# (3) Scaffold with HapHiC. 30 = haploid chromosome number (28 autosomes + X + Y).
"$HAPHIC_DIR/haphic" pipeline "$ref_file" "${output_dir}/${prefix}.filtered.bam" 30 --threads 4 --processes 2
