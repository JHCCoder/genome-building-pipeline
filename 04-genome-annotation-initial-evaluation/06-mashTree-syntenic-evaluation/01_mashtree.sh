#!/bin/bash
#SBATCH -J mashtree
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH --mem=64G
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
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_PHYLO"

# --- RUN PARAMETERS ---
REPS=100          # jackknife/bootstrap resampling replicates
NUMCPUS=8         # parallel mash sketches
# Extra arguments passed through to mash (see `mashtree --help`). Default k=21.

outdir="$SYNTENY_OUT_DIR/mashtree"
mkdir -p "$outdir"

# mashtree is a set of Perl scripts that wrap `mash` (both must be on PATH).
export PATH="$MASHTREE_DIR/bin:$(dirname "$MASH_BIN"):$PATH"

# Build jackknife + bootstrap mashtrees from the Jaccard distance between
# 21-mer sets (the tree that orders species in the synteny figure).
cat "$SYNTENY_FASTA_LIST" | xargs "$MASHTREE_DIR/bin/mashtree_jackknife.pl" \
    --reps "$REPS" --numcpus "$NUMCPUS" -- --min-depth 0 \
    > "$outdir/mashtree.jackknife.dnd"

cat "$SYNTENY_FASTA_LIST" | xargs "$MASHTREE_DIR/bin/mashtree_bootstrap.pl" \
    --reps "$REPS" --numcpus "$NUMCPUS" -- --min-depth 0 \
    > "$outdir/mashtree.bootstrap.dnd"

# Relabel leaves (raw filenames -> display species names) for the plot.
python3 "$script_dir/rename_tree.py" "$outdir/mashtree.jackknife.dnd" \
    "$SYNTENY_NAME_CONVERSIONS" "$outdir/mashtree.jackknife.renamed.dnd"
python3 "$script_dir/rename_tree.py" "$outdir/mashtree.bootstrap.dnd" \
    "$SYNTENY_NAME_CONVERSIONS" "$outdir/mashtree.bootstrap.renamed.dnd"
