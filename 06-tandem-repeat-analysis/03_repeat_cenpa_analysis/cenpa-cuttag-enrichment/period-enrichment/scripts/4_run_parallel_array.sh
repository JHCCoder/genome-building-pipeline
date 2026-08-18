#!/bin/bash
# ============================================================================
# 4_run_parallel_array.sh — Launch the MERGED-ARRAY period-enrichment pipeline,
# fully parallelized: one cluster job per bin (array 1-9) + per-sample signal.
#
#   Merge arrays      : 1b_merge_arrays.sh                  [job M]
#   Array signal      : 2_array_signal.sh        array 4   [job S]  depends M
#   Permutations      : 2d_permutation_array_per_bin.sh array 9 [job A] depends M
#   Combine+analyze   : 2f_combine_and_analyze_array.sh     [job C]  depends S+A
#
# REUSED (not re-derived): data/permutation/bedgraph/*.bedGraph (per-base CUT&Tag
# lookup tables), trf_exclusion_mask.bed, chrom_sizes, library sizes.
#
# Run from the login node:
#   sbatch scripts/4_run_parallel_array.sh   (or bash it — it only submits)
# ============================================================================

set -euo pipefail

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/scripts"

echo "[$(date)] Launching MERGED-ARRAY parallel pipeline"

M=$(sbatch "${SCRIPT_DIR}/1b_merge_arrays.sh" | awk '{print $4}')
echo "  Job M (merge arrays):      $M"

S=$(sbatch --dependency=afterok:${M} --array=0-3 "${SCRIPT_DIR}/2_array_signal.sh" | awk '{print $4}')
echo "  Job S (array signal):      $S"

A=$(sbatch --dependency=afterok:${M} --array=1-9 "${SCRIPT_DIR}/2d_permutation_array_per_bin.sh" | awk '{print $4}')
echo "  Job A (permutations 1-9):  $A"

C=$(sbatch --dependency=afterok:${S}:${A} "${SCRIPT_DIR}/2f_combine_and_analyze_array.sh" | awk '{print $4}')
echo "  Job C (combine+analyze):   $C"

echo
echo "Chain: M=${M} -> S=${S} + A=${A} -> C=${C}"
echo "Watch: watch -n 30 squeue -u jhc103"
