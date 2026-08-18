#!/bin/bash
# Submit one BISER job (01_biser_run.sh) per genome listed in BISER_GENOME_LIST.
# Run this on a login node:  bash 02_biser_submit_batch.sh

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

while read -r assembly species; do
    [ -z "$assembly" ] && continue   # skip blank lines
    echo "Submitting BISER: $assembly ($species)"
    sbatch "$script_dir/01_biser_run.sh" "$assembly" "$species"
done < "$BISER_GENOME_LIST"
