#!/bin/bash
#SBATCH -J 050525_mergeMatrix
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -t 1-00:00:00
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
conda activate "$ENV_HICEXPLORER"

genome="hifi_041425"
run="403_ear_deep_file_test"

hicMergeMatrixBins \
--matrix hicMatrix/403_ear_deep_file_test_hifi_041425_10kb.h5 --numBins 100 \
--outFileName hicMatrix/403_ear_deep_file_test_hifi_041425_10kb.100bins.h5

#(toolshed-HiCExplorer) [jhc103@tscc-1-31 hic-plot]$ seff 5109191
#Job ID: 5109191
#Cluster: tscc
#User/Group: jhc103/ren-group
#State: COMPLETED (exit code 0)
#Cores: 1
#CPU Utilized: 00:00:40
#CPU Efficiency: 76.92% of 00:00:52 core-walltime
#Job Wall-clock time: 00:00:52
#Memory Utilized: 3.33 GB
#Memory Efficiency: 33.31% of 10.00 GB
