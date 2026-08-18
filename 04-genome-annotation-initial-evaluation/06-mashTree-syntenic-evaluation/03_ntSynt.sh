#!/bin/bash
#SBATCH -J ntSynt
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -t 1-00:00:00
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
conda activate "$ENV_ASSESSMENT"

# --- RUN PARAMETERS ---
PREFIX="rodent"     # output basename for the synteny-block TSV / .fai files
DIVERGENCE=20       # approx. max % sequence divergence (ntSynt -d/--divergence)
THREADS=12

outdir="$SYNTENY_OUT_DIR/ntSynt"
mkdir -p "$outdir"
cd "$outdir"

# Compute pairwise synteny blocks across the genomes in SYNTENY_FASTA_LIST.
# Writes <PREFIX>.synteny_blocks.tsv plus one <genome>.fai per input genome
# (both are consumed by 04_ntsynt_viz.sh).
"$NTSYNT_BIN" -d "$DIVERGENCE" --fastas_list "$SYNTENY_FASTA_LIST" -p "$PREFIX" -t "$THREADS"
