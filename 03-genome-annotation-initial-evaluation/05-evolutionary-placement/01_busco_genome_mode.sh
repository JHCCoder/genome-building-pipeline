#!/bin/bash
#SBATCH -J evo_busco
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 3-00:00:00
#SBATCH --mem=160G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p platinum
#SBATCH -q hcp-csd788
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# =============================================================================
# Step 1 of the evolutionary-placement pipeline.
#
# Runs BUSCO independently on each species genome to find conserved
# single-copy orthologs, exactly as described in the methods:
#
#   busco -i <genome.fna> -m genome -l euarchontoglires_odb12 -c 16 --metaeuk
#
# One BUSCO run is produced per species (labelled by its tree_label), each in
# "$BUSCO_DIR/<label>/run_<lineage>/". Step 2 reads the resulting
# full_table.tsv files.
#
# RUN PARAMETERS
# --------------
#   BUSCO_THREADS   CPU threads handed to each BUSCO run (-c).
# =============================================================================

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

BUSCO_THREADS=16

BUSCO_DIR="$EVO_OUT_DIR/01_busco"
mkdir -p "$BUSCO_DIR"

source ~/.bashrc                 # initialize conda (adjust for your setup)
conda activate "$ENV_BUSCO"      # BUSCO + MetaEuk + Augustus

# One BUSCO run per species. The full_table.tsv and single-copy amino-acid
# sequences live under "<label>/run_<lineage>/{full_table.tsv,busco_sequences/}".
while IFS=$'\t' read -r genome label; do
    [[ -z "$genome" || "$genome" == \#* ]] && continue

    echo "==== BUSCO genome mode: $label ($genome) ===="
    busco -i "$genome" \
          -m genome \
          -l "$EVO_BUSCO_LINEAGE" \
          -c "$BUSCO_THREADS" \
          --metaeuk \
          --out_path "$BUSCO_DIR" \
          -o "$label"
done < "$EVO_SPECIES_LIST"

echo "Done. BUSCO outputs are in $BUSCO_DIR"
