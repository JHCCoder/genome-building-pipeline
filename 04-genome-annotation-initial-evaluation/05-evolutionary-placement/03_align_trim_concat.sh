#!/bin/bash
#SBATCH -J evo_align
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH --mem=120G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p platinum
#SBATCH -q hcp-csd788
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# =============================================================================
# Step 3 of the evolutionary-placement pipeline.
#
# Independently aligns (MAFFT) and trims (trimAl) the amino-acid sequences of
# each shared BUSCO gene, then concatenates all trimmed alignments into a single
# amino-acid supermatrix (step 4, via 04_concat_supermatrix.py).
#
# RUN PARAMETERS
# --------------
#   MAFFT_OPTS     extra options passed to MAFFT (default: --auto)
#   TRIMAL_MODE    trimAl trimming mode (default: -automated1)
# =============================================================================

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

MAFFT_OPTS="--auto"
TRIMAL_MODE="-automated1"

FASTA_DIR="$EVO_OUT_DIR/02_single_copy_orthologs"
ALN_DIR="$EVO_OUT_DIR/03_alignments"
GENE_LIST="$EVO_OUT_DIR/shared_singlecopy_genes.txt"
mkdir -p "$ALN_DIR"

source ~/.bashrc                      # initialize conda (adjust for your setup)
conda activate "$ENV_EVOLUTIONARY_TREE"   # MAFFT 7.526 + trimAl 1.5

# Align + trim each shared gene independently.
while IFS= read -r gene; do
    [[ -z "$gene" ]] && continue
    echo "==== aligning & trimming: $gene ===="
    mafft $MAFFT_OPTS "$FASTA_DIR/$gene.faa" > "$ALN_DIR/$gene.aln.faa"
    trimal -in "$ALN_DIR/$gene.aln.faa" \
           -out "$ALN_DIR/$gene.trim.faa" \
           $TRIMAL_MODE
done < "$GENE_LIST"

# Concatenate the trimmed alignments into the PHYLIP supermatrix.
python3 "$_script_dir/04_concat_supermatrix.py" \
    --species-list "$EVO_SPECIES_LIST" \
    --gene-list    "$GENE_LIST" \
    --trim-dir     "$ALN_DIR" \
    --out          "$EVO_OUT_DIR/04_supermatrix/supermatrix.phy"

echo "Done. Supermatrix written to $EVO_OUT_DIR/04_supermatrix/supermatrix.phy"
