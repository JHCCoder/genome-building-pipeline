#!/bin/bash
#SBATCH -J s7_coverage_query
#SBATCH -N 1 -n 1 -c 4
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

# Supplemental Figure S7 — HiFi read-coverage query for paralog vs single-copy
# gene loci. Builds the tier-classified regions (fast) then queries per-base
# HiFi coverage for each +/-5 kb locus.
#
# Usage:  sbatch 03_run_coverage.sh

set -euo pipefail

# Load shared configuration (config.sh at the repo root).
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

# Work in this script's directory (where the .py scripts and their outputs live).
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

source ~/.bashrc   # initialize conda (adjust for your setup)
# `genome-assembly` provides pandas/numpy; `02_query_coverage.py` also needs
# htslib's `tabix` via `module load` (see HTSLIB_MODULES in that script).
conda activate genome-assembly

echo "Starting region build at $(date)"
python3 01_build_regions.py

echo "Starting coverage query at $(date)"
python3 02_query_coverage.py

echo "Done at $(date)"
