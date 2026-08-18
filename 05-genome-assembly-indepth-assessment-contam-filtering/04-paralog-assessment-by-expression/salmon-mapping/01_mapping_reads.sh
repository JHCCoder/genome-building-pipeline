#!/bin/bash
#SBATCH -J salmon_mapping
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 2-12:00:00
#SBATCH --mem=120G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc
conda activate "$ENV_TRANSCRIPTOME_MAPPING"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "$SCRIPT_DIR"
SALMON_INDEX="$SCRIPT_DIR/salmon_index"

# One salmon quant per RNA-seq SRA sample
for base in SRR17216293 SRR17216294 SRR17216295 SRR17216296 SRR17216297 SRR17216298 SRR17216299 SRR17216300 SRR17216301 SRR17216302 SRR17216303 SRR17216304 SRR17216305 SRR17216306 SRR17216307 SRR17216308 SRR17216309 SRR17216310 SRR17216311 SRR17216312 SRR17216313 SRR17216314 SRR17216315 SRR17216316 SRR17216317 SRR17216318 SRR17216319 SRR17216320 SRR17216321
do
    echo "Processing $base with Salmon"
    salmon quant -i "$SALMON_INDEX" -l A \
        -1 "$MRNA_DIR/${base}_1.fastq.gz" -2 "$MRNA_DIR/${base}_2.fastq.gz" \
        -p 12 --validateMappings --gcBias -o "${base}_salmon_quant"
done

echo "All samples processed with Salmon!"
