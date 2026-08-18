#!/bin/bash
# =============================================================================
# Step 2 of the evolutionary-placement pipeline.
#
# Compares the BUSCO full_table.tsv outputs across all species, retains the
# single-copy genes classified "Complete" in every genome, and writes one
# multi-species FASTA per retained gene (logic in
# 02_collect_shared_singlecopy_orthologs.py).
#
# Lightweight (parsing + small file copies) — run directly, no Slurm needed.
# Must run after 01_busco_genome_mode.sh finishes.
# =============================================================================

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

python3 "$_script_dir/02_collect_shared_singlecopy_orthologs.py" \
    --species-list "$EVO_SPECIES_LIST" \
    --busco-dir    "$EVO_OUT_DIR/01_busco" \
    --lineage      "$EVO_BUSCO_LINEAGE" \
    --out-dir      "$EVO_OUT_DIR/02_single_copy_orthologs" \
    --gene-list    "$EVO_OUT_DIR/shared_singlecopy_genes.txt"
