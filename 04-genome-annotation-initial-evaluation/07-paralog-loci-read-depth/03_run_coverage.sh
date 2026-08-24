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

# Locate this script's directory.
#
# Under Slurm, BASH_SOURCE[0] points at the *spool copy* of the script
# (/cm/local/.../slurm_script), not at the file you submitted, so resolving the
# repo root from it walks up to "/" and fails. SLURM_SUBMIT_DIR is the directory
# sbatch was invoked from, so prefer it when running as a batch job and fall
# back to BASH_SOURCE for a plain `bash 03_run_coverage.sh`.
if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    _script_dir="$SLURM_SUBMIT_DIR"
else
    _script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
fi

# Load shared configuration (config.sh at the repo root).
_repo_root="$_script_dir"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
if [[ ! -f "$_repo_root/config.sh" ]]; then
    echo "ERROR: could not locate config.sh above $_script_dir" >&2
    echo "       (submit with 'cd <this script's directory> && sbatch 03_run_coverage.sh')" >&2
    exit 1
fi
source "$_repo_root/config.sh"
unset _repo_root

# Work in this script's directory (where the .py scripts and their outputs live).
cd "$_script_dir" || exit 1
unset _script_dir

source ~/.bashrc   # initialize conda (adjust for your setup)
# `genome-assembly` provides pandas/numpy; `02_query_coverage.py` also needs
# htslib's `tabix` via `module load` (see HTSLIB_MODULES in that script).
#
# `set -u` must be lifted across conda activation: some packages ship
# activate.d hooks that reference unset variables (here openjdk_activate.sh
# reads $target_platform), which aborts the job under `set -euo pipefail`.
set +u
conda activate genome-assembly
set -u

echo "Starting region build at $(date)"
python3 01_build_regions.py

echo "Starting coverage query at $(date)"
python3 02_query_coverage.py

echo "Done at $(date)"
