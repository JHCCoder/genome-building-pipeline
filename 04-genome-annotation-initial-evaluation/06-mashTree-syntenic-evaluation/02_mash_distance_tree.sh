#!/bin/bash
#SBATCH -J mash_dist_tree
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 8:00:00
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
conda activate "$ENV_PHYLO"   # provides fastme (mash is a standalone binary)

# --- RUN PARAMETERS ---
SKETCH_SIZE=8000   # hash count per genome
THREADS=8

outdir="$SYNTENY_OUT_DIR/mash-distance"
mkdir -p "$outdir/sketches"
cd "$outdir"

# 1. Sketch each genome (default k=21; Jaccard distance is computed over 21-mers).
while IFS= read -r f; do
  "$MASH_BIN" sketch "$f" -s "$SKETCH_SIZE" -p "$THREADS"
  mv "$f.msh" sketches/
done < "$SYNTENY_FASTA_LIST"

# 2. All-vs-all Jaccard distances (lower-triangular matrix, leading taxon count).
"$MASH_BIN" triangle sketches/*.msh > all_distances.tab

# 3. Strip directory prefixes from the row labels.
sed 's|/.*/||g' all_distances.tab > all_distances_cleaned.tab

# 4. Lower-triangular -> square PHYLIP matrix (FastME input).
python3 "$script_dir/convert_to_square.py" all_distances_cleaned.tab all_distances_cleaned.phy

# 5. Distance tree (BIONJ) from the Jaccard distance matrix.
fastme -i all_distances_cleaned.phy -o all_distances_cleaned.nwk -D 101

# 6. Relabel leaves (raw filenames -> display species names).
python3 "$script_dir/rename_tree.py" all_distances_cleaned.nwk \
    "$SYNTENY_NAME_CONVERSIONS" all_distances_cleaned_renamed.nwk
