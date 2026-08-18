#!/bin/bash
#SBATCH -J purge_dup_align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 2-00:00:00
#SBATCH --mem=48G
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
conda activate "$ENV_ASSESSMENT"

pri_asm="$CHR_ASSIGNED_ASSEMBLY"
mkdir -p "$PURGE_WORK_DIR"
cd "$PURGE_WORK_DIR"

# Align HiFi reads to the assembly (minimap2 asm20), then compute coverage stats
for i in $PURGE_PB_LIST; do
    file_name=$(basename "$i")
    new_file_name="${file_name%.fastq.gz}.paf.gz"
    minimap2 -xasm20 -I -t 8 "$pri_asm" "$i" | gzip -c - > "$new_file_name"
done

"$PURGE_DUPS_BIN/pbcstat" *.paf.gz
"$PURGE_DUPS_BIN/calcuts" PB.stat > cutoffs 2>calcults.log

"$PURGE_DUPS_BIN/split_fa" "$pri_asm" > "$pri_asm.split"
minimap2 -xasm5 -DP "$pri_asm.split" "$pri_asm.split" | gzip -c - > "$pri_asm.split.self.paf.gz"
