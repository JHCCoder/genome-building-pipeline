#!/bin/bash
# ============================================================================
# 9_run_chain.sh -- orchestrate the CENP-A domain-centered repeat composition
# pipeline with Slurm dependencies.
#
# Stages and dependencies:
#   1_nh_weighted    (array 150/151/152)   [running or done]
#   3a_9bin_prep                           [done]
#   3_covariates    (mapp/GC/repdensity)  depends 3a
#   2_domains       (R, interactive)      depends 1   -> domains
#   4_composition   (R, interactive)      depends 2 + 3a
#   5_matched_null  (Python, slurm)       depends 3 + 2
#   5b_null_analysis(R, interactive)      depends 5
#   6a_extract_seqs                       [submitted]
#   6b_kmc_count    depends 6a
#   6c_probe_setup  depends 6b
#   6d_count_probes depends 6c
#   6e_kmer_analyze (R, interactive)      depends 6d
#   7_signal        (R, interactive)      depends 1 + 2
#   8_plots         (R, interactive)      depends 4 + 5b + 6e + 7
#
# Usage: bash scripts/9_run_chain.sh   (must have stage1 array job ID set)
# ============================================================================

set -euo
cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh

STEP1=${STEP1_JOB:-}
echo "STEP1 job id: ${STEP1:-<unset>}"

# If step1 not given, find it in squeue / sacct
if [[ -z "$STEP1" ]]; then
    STEP1=$(squeue -u jhc103 -h -o "%i" -n 0807_nh_weighted 2>/dev/null | head -1 || true)
    echo "auto-detected STEP1=$STEP1"
fi

DEP1=""
if [[ -n "$STEP1" ]]; then DEP1="--dependency=afterok:${STEP1}"; fi

# ---- 6b KMC count (depends 6a) ----
JOB6A=$(squeue -u jhc103 -h -o "%i" -n 0807_kmer_seqs 2>/dev/null | head -1 || true)
echo "6a job: ${JOB6A:-<done-or-running>}"
if [[ -n "$JOB6A" ]]; then
    JOB6B=$(sbatch --dependency=afterok:${JOB6A} scripts/6b_kmc_count.sh | awk '{print $4}')
else
    JOB6B=$(sbatch scripts/6b_kmc_count.sh | awk '{print $4}')
fi
echo "6b -> $JOB6B"
JOB6C=$(sbatch --dependency=afterok:${JOB6B} scripts/6c_probe_setup.sh | awk '{print $4}')
echo "6c -> $JOB6C"
JOB6D=$(sbatch --dependency=afterok:${JOB6C} scripts/6d_count_probes.sh | awk '{print $4}')
echo "6d -> $JOB6D"

# ---- 5 matched null (depends 3 covariates + 2 domains) ----
# 3_covariates and 3a must finish first
JOB3=$(squeue -u jhc103 -h -o "%i" -n 0807_covariates 2>/dev/null | head -1 || true)
echo "3 job: ${JOB3:-<done-or-running>}"
JOB3A=$(squeue -u jhc103 -h -o "%i" -n 0807_9bin_prep 2>/dev/null | head -1 || true)
echo "3a job: ${JOB3A:-<done>}"
DEP5=""
if [[ -n "$JOB3" ]]; then DEP5="--dependency=afterok:${JOB3}"; fi
JOB5=$(sbatch $DEP5 scripts/5_matched_null.sh | awk '{print $4}')
echo "5 -> $JOB5"

echo "=== Interactive stages (run after deps): 2_domains.R, 4_composition_9bin.R, 5b_null_analysis.R, 6e_kmer_analyze.R, 7_signal.R, 8_plots.R ==="
echo "Suggested chain (after stage 1 array + 6a finish):"
echo "  conda activate r-visualizations"
echo "  Rscript scripts/2_domains.R <workdir> weighted"
echo "  Rscript scripts/2_domains.R <workdir> k1"
echo "  Rscript scripts/4_composition_9bin.R <workdir>"
echo "  Rscript scripts/5b_null_analysis.R <workdir>"
echo "  Rscript scripts/6e_kmer_analyze.R <workdir>"
echo "  Rscript scripts/7_signal.R <workdir>"
echo "  Rscript scripts/8_plots.R <workdir>"
