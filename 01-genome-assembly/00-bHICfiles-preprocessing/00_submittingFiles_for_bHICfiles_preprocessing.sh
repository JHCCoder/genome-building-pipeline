#!/bin/bash
# Submit one 00_bHICfiles_preprocessing.sh job per sample listed in the Hi-C file
# list. Run directly (e.g. `bash 00_submittingFiles_for_bHICfiles_preprocessing.sh`).

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

echo "Starting job submission from $BHIC_FILE_LIST ..."
count=0

while IFS=$'\t' read -r read_file1 read_file2 sampleid; do
  count=$((count + 1))
  echo "Processing job $count:"
  echo "  R1: $read_file1"
  echo "  R2: $read_file2"
  echo "  Sample: $sampleid"

  job_id=$(sbatch "$SCRIPT_DIR/00_bHICfiles_preprocessing.sh" "$read_file1" "$read_file2" "$sampleid")
  echo "  Submitted job: $job_id"
  echo "---"
done < "$BHIC_FILE_LIST"

echo "Total jobs submitted: $count"
