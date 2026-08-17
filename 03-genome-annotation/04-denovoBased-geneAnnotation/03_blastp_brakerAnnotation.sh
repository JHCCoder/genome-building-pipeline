#!/bin/bash
#SBATCH -J blastp_braker
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH --mem=2G
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

# NOTE: this is a benchmark run on a 100-entry subset (braker_121624_100_entries.aa);
# the full functional-annotation run used a parallelized blastp script.
echo "running blastp"
export BLASTDB="$BLAST_DB_DIR"
"$BLASTP_BIN" -query braker_121624_100_entries.aa -db UniProtVertebrate_DB \
  -out "$SCRATCH_DIR/output-blastp-test/braker_121624_100_entries_vertebrateUniprot_database_thread1.tsv" \
  -evalue 1e-5 -num_threads 1 -max_target_seqs 8 \
  -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -seg yes
