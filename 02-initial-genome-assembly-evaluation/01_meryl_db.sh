#!/bin/bash
#SBATCH -J meryl_kmer_db
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -t 12:00:00
#SBATCH --mem=128G
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
conda activate "$ENV_MERQURY"

mkdir -p "$MERQURY_OUT_DIR"
cd "$MERQURY_OUT_DIR"

# Build a k=21 k-mer count database from Illumina WGS reads (Merqury prerequisite).
# NOTE: 02_merqury.sh consumes the db named by $MERYL_DB (collaborator_degu_WGS.meryl);
# rename or symlink this output to match.
meryl k=21 count "$WGS_READS" output male403_WGS_collaborator.meryl
