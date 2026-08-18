#!/bin/bash
#SBATCH -J 071626_verify_diag
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 2:00:00
#SBATCH --mem=16G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com

source /tscc/nfs/home/jhc103/.bashrc
conda activate r-visualizations

Rscript verify_diagnostics.R
