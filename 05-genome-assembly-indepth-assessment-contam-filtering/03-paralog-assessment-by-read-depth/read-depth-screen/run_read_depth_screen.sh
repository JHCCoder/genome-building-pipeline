#!/bin/bash
#SBATCH -J read_depth_screen
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH -t 2-00:00:00
#SBATCH --mem=256G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc
# Temporarily disable -u so conda's openjdk_activate.sh doesn't crash
set +u
conda activate "$ENV_ASSESSMENT"
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

ASSEMBLY="$READ_DEPTH_ASSEMBLY"
HIFI_READS="$READ_DEPTH_HIFI_DIR"
PARALOG_FAMILIES="$SCRIPT_DIR/../paralog_families.tsv"   # from build_paralog_families.py
REPEAT_BED="$PARALOG_REPEAT_BED"
GENOME_LENGTH_FILE="$PARALOG_GENOME_LENGTHS"

BAM_PREFIX="hifi_male_map-hifi_contamFiltered_chrAssigned"
OUT_DIR="${SCRIPT_DIR}"
mkdir -p "${OUT_DIR}"

# Step 1: map HiFi reads with minimap2 map-hifi → sort → index
echo "STEP 1: mapping HiFi reads with minimap2 map-hifi"
BAM_FILE="${OUT_DIR}/${BAM_PREFIX}.bam"
if [ -f "${BAM_FILE}" ] && [ -f "${BAM_FILE}.bai" ]; then
    echo "BAM and index already exist: ${BAM_FILE} — skipping mapping."
else
    READ_FILES=()
    for f in "${HIFI_READS}"/*.fastq.gz; do
        [ -f "$f" ] && READ_FILES+=("$f")
    done
    echo "Found ${#READ_FILES[@]} HiFi read files."
    # -F 0x900 excludes supplementary (0x800) and secondary (0x100) alignments
    minimap2 -ax map-hifi -t 36 "${ASSEMBLY}" "${READ_FILES[@]}" \
        | samtools view -@ 8 -F 0x900 -Sb \
        | samtools sort -@ 8 -o "${BAM_FILE}" \
        && samtools index -@ 8 "${BAM_FILE}"
fi

# Step 2: coverage analysis / haplotig screening
echo "STEP 2: coverage analysis and haplotig screening"
python3 "${SCRIPT_DIR}/read_depth_screen.py" \
    --bam "${BAM_FILE}" \
    --fasta "${ASSEMBLY}" \
    --paralog-families "${PARALOG_FAMILIES}" \
    --genome-lengths "${GENOME_LENGTH_FILE}" \
    --out-dir "${OUT_DIR}" \
    --repeat-bed "${REPEAT_BED}" \
    --threads 24 \
    --flank-bp 25000

echo "PIPELINE COMPLETE — output: ${OUT_DIR}"
