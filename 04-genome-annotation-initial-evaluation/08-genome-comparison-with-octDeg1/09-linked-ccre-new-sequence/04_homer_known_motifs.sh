#!/bin/bash
#SBATCH -J 101326_homer_linkedPeaks_notShared
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --mem=48G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user=you@example.com   # TODO: your email
#SBATCH --no-requeue

# Step 4: HOMER known-motif enrichment in the linked peaks NOT shared with OctDeg1.0
# (step 3 classes new_sequence + partial + resolved_repeat) against the linked peaks shared
# with OctDeg1.0 (class shared) as background. HOMER matches GC content of the background.
#
#   -size 200   200 bp centred on each peak (Cell Ranger ARC reports no summits; peaks are
#               ~870 bp, so the full width would dilute motif signal)
#   -nomotif    known motifs only (HOMER vertebrate library)
#   no -mask    the OctDeg2.0 FASTA is not soft-masked, and repeat-derived motifs are part of
#               the question; interpret transposon-associated motifs accordingly
#
# Outputs: $FIG3D_CCRE_DIR/04_homer_notShared_vs_shared/ (knownResults.txt + html)
#          results/04_homer_knownResults_notShared_vs_shared.txt (copy)
set -euo pipefail

if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    _script_dir="$SLURM_SUBMIT_DIR"
else
    _script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
fi
_repo_root="$_script_dir"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root
export PATH="$HOMER_BIN_DIR:$PATH"

OUT="$FIG3D_CCRE_DIR"
HDIR="$OUT/04_homer_notShared_vs_shared"
[[ -e "$HDIR" ]] && { echo "ERROR: $HDIR exists; move it aside before re-running" >&2; exit 1; }
for c in shared new_sequence partial resolved_repeat; do
  [[ -s "$OUT/03_linked_peaks_${c}.bed" ]] || { echo "ERROR: run 03_classify_new_sequence.py first" >&2; exit 1; }
done

# HOMER peak BED: chr start end unique_id score strand
TARGET="$OUT/04_target_notShared.bed"
BG="$OUT/04_background_shared.bed"
cat "$OUT"/03_linked_peaks_{new_sequence,partial,resolved_repeat}.bed \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".","+"}' > "$TARGET"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,".","+"}' "$OUT/03_linked_peaks_shared.bed" > "$BG"
echo "target peaks: $(wc -l < "$TARGET")   background peaks: $(wc -l < "$BG")"

findMotifsGenome.pl "$TARGET" "$MITO_ASSEMBLY" "$HDIR" \
  -bg "$BG" -size 200 -nomotif -p 8 -preparsedDir "$OUT/04_homer_preparsed"

mkdir -p "$_script_dir/results"
cp "$HDIR/knownResults.txt" "$_script_dir/results/04_homer_knownResults_notShared_vs_shared.txt"
echo "homer finished: $(date)"
