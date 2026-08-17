#!/bin/bash
# ============================================================================
# 4_run_parallel.sh — Launch the re-binned (9-bin) period-enrichment pipeline,
# fully parallelized: one cluster job per bin (array 1-9), plus shared setup.
#
#   Step 1 (rebin)      : 1_prepare_period_bins.sh          [job P]
#   Setup (bedGraph)    : 2e_prepare_bedgraphs.sh  array 4  [job BG]
#   Permutations        : 2d_permutation_per_bin.sh array 9 [job PBIN]  depends P+BG
#   Combine+analyze     : 2f_combine_and_analyze.sh         [job CA]    depends PBIN
#
# Run from the login node:
#   sbatch scripts/4_run_parallel.sh   (or bash it — it only submits)
# ============================================================================

set -euo pipefail

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/scripts"

echo "[$(date)] Launching re-binned 9-bin parallel pipeline"

# Step 1: re-bin TRF intervals
P=$(sbatch "${SCRIPT_DIR}/1_prepare_period_bins.sh" | awk '{print $4}')
echo "  Job P (rebin):        $P"

# Setup: shared bedGraphs (1 job per sample, 4 tasks)
BG=$(sbatch "${SCRIPT_DIR}/2e_prepare_bedgraphs.sh" | awk '{print $4}')
echo "  Job BG (bedGraphs):   $BG"

# Per-bin permutations (9 parallel jobs)
PBIN=$(sbatch --dependency=afterok:${P}:${BG} --array=1-9 \
      "${SCRIPT_DIR}/2d_permutation_per_bin.sh" | awk '{print $4}')
echo "  Job PBIN (array 1-9): $PBIN"

# Combine + re-tag + analysis + violin (after all bins finish)
CA=$(sbatch --dependency=afterok:${PBIN} "${SCRIPT_DIR}/2f_combine_and_analyze.sh" | awk '{print $4}')
echo "  Job CA (combine+analysis): $CA"

echo
echo "Chain: P=${P} BG=${BG} -> PBIN=${PBIN} -> CA=${CA}"
echo "Watch: watch -n 30 squeue -u jhc103"
