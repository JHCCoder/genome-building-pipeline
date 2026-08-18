#!/bin/bash
#SBATCH -J ntsynt_viz
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --mem=32G
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
conda activate "$ENV_ASSESSMENT"   # provides snakemake (ntSynt-viz backend)

# --- RUN PARAMETERS ---
RIBBON_ADJUST=0.14
SCALE=1e9

blocks="$SYNTENY_OUT_DIR/ntSynt/rodent.synteny_blocks.tsv"
fais="$SYNTENY_OUT_DIR/ntSynt"/*.fai

# Ribbon plot of OctDeg2.0 (target genome) against the other species, with the
# species ordered along the y-axis by the mashtree from 01_mashtree.sh. Two
# variants: the bootstrap and jackknife mashtrees (same topology, different
# support values).
"$NTSYNT_VIZ_BIN" \
  --blocks "$blocks" \
  --fais $fais \
  --normalize \
  --prefix rodent_ribbon_plot_mashTree_bootStrap \
  --ribbon_adjust "$RIBBON_ADJUST" \
  --scale "$SCALE" \
  --name_conversion "$SYNTENY_NAME_CONVERSIONS" \
  --target-genome "$SYNTENY_TARGET_GENOME" \
  --tree "$SYNTENY_OUT_DIR/mashtree/mashtree.bootstrap.renamed.dnd"

"$NTSYNT_VIZ_BIN" \
  --blocks "$blocks" \
  --fais $fais \
  --normalize \
  --prefix rodent_ribbon_plot_mashTree_jackknife \
  --ribbon_adjust "$RIBBON_ADJUST" \
  --scale "$SCALE" \
  --name_conversion "$SYNTENY_NAME_CONVERSIONS" \
  --target-genome "$SYNTENY_TARGET_GENOME" \
  --tree "$SYNTENY_OUT_DIR/mashtree/mashtree.jackknife.renamed.dnd"
