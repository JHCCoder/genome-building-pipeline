#!/bin/bash
#SBATCH -J 082225_hicBuildMatrix
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=178G
#SBATCH -t 7-00:00:00
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p hotel
#SBATCH -q hotel
#SBATCH -A htl195
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
run="403_ear_deep_file"

hicBuildMatrix --samFiles ${run}_${genome}_1.sorted.bam ${run}_${genome}_2.sorted.bam --binSize 10000 --restrictionSequence GATC --danglingSequence GATC --restrictionCutFile $HIC_RESTRICTION_BED --outBam ${genome}_ref.bam --outFileName hicMatrix-longRead/${run}_${genome}_10kb.h5 --QCfolder hicMatrix-longRead/${run}_${genome}_10kb_QC --threads 4 --inputBufferSize 400000


#(toolshed-HiCExplorer) [jhc103@tscc-1-31 hic-plot]$ seff 5107396
#Job ID: 5107396
#Cluster: tscc
#User/Group: jhc103/ren-group
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 8
#CPU Utilized: 06:53:53
#CPU Efficiency: 49.94% of 13:48:48 core-walltime
#Job Wall-clock time: 01:43:36
#Memory Utilized: 198.49 GB
#Memory Efficiency: 103.38% of 192.00 GB
