#!/bin/bash
#SBATCH -J repeatModeler
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 33
#SBATCH -t 7-00:00:00
#SBATCH --mem=660G
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
conda activate "$ENV_REPEATMODELER"

# Build a de-novo repeat library (RepeatModeler) from the scaffolded assembly,
# then mask with RepeatMasker (see 02_repeatMasker.sh).
work_dir="$REPEAT_OUT_DIR/hifiasm-121624-haphic"
mkdir -p "$work_dir"
cd "$work_dir"

# Step 1: build a BLAST database of the assembly
BuildDatabase -name odegus "$REPEAT_MODELER_ASM"

# Step 2: build the de-novo repeat library (~17-18 h)
RepeatModeler -database odegus -threads 32 -LTRStruct -engine ncbi > out.log
