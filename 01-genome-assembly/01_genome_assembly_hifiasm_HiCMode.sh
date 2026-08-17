#!/bin/bash
#SBATCH -J hifiasm_hic_mode
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH -t 4-00:00:00
#SBATCH --mem=256G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p platinum
#SBATCH -q hcp-csd788
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

mkdir -p "$HIFIASM_OUT_DIR"

# Assemble HiFi reads with hifiasm (Hi-C + HiFi mode; --primary keeps one haplotype).
# For a hybrid ONT+HiFi assembly, append the ONT reads to the input list as well.
"$HIFIASM_BIN" -o "$HIFIASM_OUT_DIR/$HIFIASM_ASM_NAME" \
  --primary -t36 -l3 -k21 \
  --h1 "$HIC_RAW_R1" \
  --h2 "$HIC_RAW_R2" \
  "$HIFI_READS"
