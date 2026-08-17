#!/bin/bash
#SBATCH -J salmon_index
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 4:00:00
#SBATCH --mem=64G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue

# Selective-alignment salmon index (combine-lab alevin tutorial)

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc
conda activate "$ENV_TRANSCRIPTOME_MAPPING"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "$SCRIPT_DIR"

# Build the transcriptome + decoys from the annotation GFF + assembly, then index
gffread "$LIFTOFF_REF_GFF" -g "$MITO_ASSEMBLY" -w all_transcripts.fa
grep "^>" "$MITO_ASSEMBLY" | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat all_transcripts.fa "$MITO_ASSEMBLY" > transcriptome_genome.fa

salmon index -t transcriptome_genome.fa -d decoys.txt -p 12 -i salmon_index --gencode
