#!/bin/bash
#SBATCH -J 101226_cellrangerArc_mkref_octDeg2_liftoffOnly
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 6:00:00
#SBATCH --mem=64G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue

# Build the Cell Ranger ARC reference for OctDeg2.0 with Liftoff-only gene models
# (OctDeg1.0 RefSeq lifted onto OctDeg2.0, matched to the OctDeg1.0 reference by
# 05_build_liftoff_only_gtf.py). Panel D then compares the assemblies on matched gene
# models, without BRAKER genes or peaks2utr UTRs.
#
# Settings mirror the OctDeg1.0 reference ($FIG3D_OCTDEG1_CR_REF): no TF motif file, and
# no non_nuclear_contigs (OctDeg1.0 contains the mitogenome NC_020661.1 but did not flag
# it). Note that the earlier OctDeg2 references DID flag NC_020661.1.

set -euo pipefail

# Under sbatch, BASH_SOURCE is Slurm's spool copy of this script, so walking up from it
# never reaches the repo. SLURM_SUBMIT_DIR is where sbatch was invoked; prefer it in a
# batch job and fall back to BASH_SOURCE for a plain `bash 06_...sh`.
if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    _script_dir="$SLURM_SUBMIT_DIR"
else
    _script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
fi

# Load shared configuration (config.sh at the repo root)
_repo_root="$_script_dir"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

GENOME="$(basename "$FIG3D_OCTDEG2_CR_REF")"          # OctDeg2_liftoffOnly
if [[ -e "$FIG3D_OCTDEG2_CR_REF" ]]; then
  echo "ERROR: $FIG3D_OCTDEG2_CR_REF already exists; move it aside before rebuilding" >&2
  exit 1
fi
[[ -s "$FIG3D_LIFTOFF_GTF" ]] || { echo "ERROR: run 05_build_liftoff_only_gtf.py first" >&2; exit 1; }

mkdir -p "$FIG3D_CR_DIR"
cd "$FIG3D_CR_DIR"            # mkref writes ./$GENOME

CONFIG="$FIG3D_CR_DIR/${GENOME}.mkref.config"
cat > "$CONFIG" <<EOF
{
    organism: "Octodon_degus"
    genome: ["$GENOME"]
    input_fasta: ["$MITO_ASSEMBLY"]
    input_gtf: ["$FIG3D_LIFTOFF_GTF"]
}
EOF
cat "$CONFIG"

echo "cellranger-arc: $("$CELLRANGER_ARC_BIN" --version 2>&1 | head -1)"
"$CELLRANGER_ARC_BIN" mkref --config="$CONFIG" --nthreads=16 --memgb=56

echo "mkref finished: $(date)"
