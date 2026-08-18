#!/bin/bash
# Build the STAR genome index for transcriptome mapping (splice-aware).

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc
conda activate "$ENV_TRANSCRIPTOME_MAPPING"

STAR --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX" \
     --genomeFastaFiles "$MITO_ASSEMBLY" \
     --sjdbGTFfile cct7_complete_annotations.gff \
     --genomeSAindexNbases 12
