#!/bin/bash
# ============================================================
# Submit HiCAT jobs for ALL chromosomes with tiered resources
# ============================================================
#
# Resource tiers based on chr4:125-135mb (10Mb) baseline:
#   that run used 159 GB memory, 18h wall time on 8 cores
#
# Tiers scale by chromosome size, using 32 cores for all.
#
# Usage: bash submit_HiCAT_all.sh
#   Set DRY_RUN=1 to print commands without submitting.
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUNNER="${SCRIPT_DIR}/run_HiCAT_chr.sh"

DRY_RUN="${DRY_RUN:-0}"

# --------------------------------------------------
# Tier definitions: chr_list | time | memory
# --------------------------------------------------
# Tier 1: Largest chromosomes (~150-200 Mb)
TIER1_CHRS=(chrX chr1 chr2)
TIER1_TIME="7-00:00:00"
TIER1_MEM="300G"

# Tier 2: Large chromosomes (~125-150 Mb)
TIER2_CHRS=(chr3 chr4 chr5 chr6 chr7 chr8 chr9)
TIER2_TIME="5-00:00:00"
TIER2_MEM="250G"

# Tier 3: Medium chromosomes (~100-125 Mb)
TIER3_CHRS=(chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17)
TIER3_TIME="3-00:00:00"
TIER3_MEM="200G"

# Tier 4: Small chromosomes (<100 Mb)
TIER4_CHRS=(chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chrY)
TIER4_TIME="2-00:00:00"
TIER4_MEM="150G"

submit_tier() {
    local -n chrs=$1
    local time=$2
    local mem=$3

    for chr in "${chrs[@]}"; do
        local jobname="HiCAT_${chr}"
        local cmd="sbatch \
            --job-name=${jobname} \
            --time=${time} \
            --mem=${mem} \
            ${RUNNER} ${chr}"

        echo "[submit] ${cmd}"
        if [[ "${DRY_RUN}" != "1" ]]; then
            eval "${cmd}"
        fi
    done
}

echo "=== HiCAT genome-wide submission ==="
echo "Runner:  ${RUNNER}"
echo "Mode:    $([ "${DRY_RUN}" == "1" ] && echo 'DRY RUN' || echo 'LIVE')"
echo ""

submit_tier TIER1_CHRS "${TIER1_TIME}" "${TIER1_MEM}"
submit_tier TIER2_CHRS "${TIER2_TIME}" "${TIER2_MEM}"
submit_tier TIER3_CHRS "${TIER3_TIME}" "${TIER3_MEM}"
submit_tier TIER4_CHRS "${TIER4_TIME}" "${TIER4_MEM}"

echo ""
echo "Done. Submitted 30 chromosomes total."
echo "Outputs will go to: ${SCRIPT_DIR}/HiCAT_genome/<chr>/"
echo ""
echo "Monitor with: squeue -u jhc103 | grep HiCAT"
echo "Or:          seff <jobid>"
