#!/bin/bash
#SBATCH -J inspector
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=10G
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
conda activate "$ENV_INSPECTOR"

# --- RUN PARAMETERS ---
THREADS=8
DATATYPE="hifi"          # PacBio HiFi reads (--datatype clr|hifi|nanopore)
OUT_SUBDIR="hifiasm-041425-scaffold"

mkdir -p "$INSPECTOR_OUT_DIR/$OUT_SUBDIR"

# Reference-free evaluation: maps HiFi reads to the assembly and reports structural
# errors (expansion/collapse/haplotype-switch/inversion), small-scale errors and QV.
# NOTE: -r is left unquoted so the wildcard in $INSPECTOR_HIFI_READS expands.
"$INSPECTOR_BIN" -c "$INSPECTOR_ASM" -r $INSPECTOR_HIFI_READS \
    -o "$INSPECTOR_OUT_DIR/$OUT_SUBDIR" --datatype "$DATATYPE" -t "$THREADS"

# Optional — reference-based mode (adds contig-to-reference coverage + errors).
# "$INSPECTOR_BIN" -c "$INSPECTOR_ASM" -r $INSPECTOR_HIFI_READS --ref "$INSPECTOR_REF" \
#     -o "$INSPECTOR_OUT_DIR/hifiasm-041425-scaffold-reference" --datatype "$DATATYPE" -t "$THREADS"
