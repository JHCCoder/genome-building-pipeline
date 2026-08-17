#!/bin/bash
#SBATCH -J submit_hisat2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 02:00:00
#SBATCH --mem=32G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_GENOME_ANNOTATION"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

# Build the HISAT2 index from the masked assembly (if not already present)
if [[ ! -f "${HISAT2_INDEX_NAME}.1.ht2" ]]; then
    echo "HISAT2 index not found. Building index..."
    hisat2-build -p 12 --verbose "$MASKED_ASSEMBLY" "$HISAT2_INDEX_NAME"
    if [[ $? -ne 0 ]]; then
        echo "Error: HISAT2 index build failed."
        exit 1
    fi
else
    echo "HISAT2 index already exists. Skipping index build."
fi

# Submit one array job per SRA accession (RNA-seq)
SRR_list=(
    SRR17216301 SRR17216302 SRR17216303 SRR17216304 SRR17216305
    SRR17216306 SRR17216307 SRR17216308 SRR17216309 SRR17216310
    SRR17216311 SRR17216312 SRR17216313 SRR17216314 SRR17216315
    SRR17216293 SRR17216294 SRR17216295 SRR17216296 SRR17216297
    SRR17216320 SRR17216319 SRR17216318 SRR17216317 SRR17216316
    SRR17216299 SRR17216298 SRR17216321 SRR17216300
)

sbatch --array=0-$((${#SRR_list[@]} - 1)) "$SCRIPT_DIR/01_align_transcriptome.sh" "${SRR_list[@]}"
