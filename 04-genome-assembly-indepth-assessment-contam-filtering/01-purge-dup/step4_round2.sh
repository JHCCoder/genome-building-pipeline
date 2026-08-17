#!/bin/bash
#SBATCH -J purge_dup_round2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 6:00:00
#SBATCH --mem=128G
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
hap_asm="$PURGE_HAP_ASM"
cd "$PURGE_WORK_DIR"

outdir2="$PURGE_WORK_DIR/round2"
mkdir -p "$outdir2"

# Combine hap.fa (from step3 get_seqs) with the haplotig assembly, then rerun
cat hap.fa "$hap_asm" > "$outdir2/round2.fa"
cd "$outdir2"

for i in $PURGE_PB_LIST; do
    filename=$(basename "$i")
    new_file_name="${filename%.fastq.gz}"
    minimap2 -xasm20 -t 8 round2.fa "$i" | gzip -c - > "$new_file_name.paf.gz"
done

"$PURGE_DUPS_BIN/pbcstat" *.paf.gz
"$PURGE_DUPS_BIN/calcuts" PB.stat > cutoffs 2>calcults.log

"$PURGE_DUPS_BIN/split_fa" round2.fa > round2.split
minimap2 -xasm5 -DP round2.split round2.split | gzip -c - > round2.split.self.paf.gz

"$PURGE_DUPS_BIN/purge_dups" -2 -T cutoffs -c PB.base.cov round2.split.self.paf.gz > dups.bed 2> purge_dups.log

"$PURGE_DUPS_BIN/get_seqs" -e dups.bed round2.fa
