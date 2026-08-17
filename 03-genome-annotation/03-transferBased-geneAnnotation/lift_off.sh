#!/bin/bash
#SBATCH -J liftoff_transfer
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 01:00:00
#SBATCH --mem=80G
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
conda activate "$ENV_LIFTOFF"

assembly="verkko-hap1-transferred-from-hifiasm-041425"
query_fa="$LIFTOFF_QUERY"
ref_fa="$LIFTOFF_REF"
ref_gff="$LIFTOFF_REF_GFF"
output_dir="$OUTPUT_DIR/outputs-from-liftoff/$assembly"
temp_dir="$SCRATCH_DIR/temp-dir-fast/liftoff-intermediate-files/$assembly"
mkdir -p "$output_dir" "$temp_dir"

# Transfer the reference annotation onto the de-novo assembly
liftoff "$query_fa" "$ref_fa" -g "$ref_gff" \
  -o "$output_dir/$assembly.gff" \
  -u "$output_dir/unmapped_features.txt" \
  -dir "$temp_dir" -p 4 > "liftoff_${assembly}_run.out"
