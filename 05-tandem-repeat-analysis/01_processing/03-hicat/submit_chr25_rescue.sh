#!/bin/bash
# ============================================================
# Submit chr25 rescue jobs (2-way split sub-sub-chunks).
#
# Each incomplete sub-chunk split into a/b halves:
#   chr25_part01_sub06a, chr25_part01_sub06b
#
# Conservative settings:
#   300G memory, 2 threads, 2 days, GPU nodes excluded
#
# Usage: bash submit_chr25_rescue.sh [--dry-run]
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GENOME_DIR="${SCRIPT_DIR}/HiCAT_genome/chr25"

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

# ── GPU nodes to exclude (all nodes with "gpu" in name) ───────────
GPU_NODES=$(scontrol show nodes 2>/dev/null \
    | grep -i "NodeName=" \
    | grep -i "gpu" \
    | grep -oP 'NodeName=\K\S+' \
    | sort -u \
    | paste -sd, -)

if [[ -z "${GPU_NODES}" ]]; then
    echo "WARNING: Could not determine GPU node list. Not excluding any nodes."
    EXCLUDE_ARG=""
else
    GPU_COUNT=$(echo "${GPU_NODES}" | tr ',' '\n' | wc -l)
    echo "Excluding ${GPU_COUNT} GPU nodes"
    EXCLUDE_ARG="--exclude=${GPU_NODES}"
fi

# ── Resources ─────────────────────────────────────────────────────
MEM="300G"
CPUS=2

# ── Discover sub-sub-chunk fastas ─────────────────────────────────
echo "=== Discovering rescue sub-sub-chunk fastas ==="
RESCUE_CHUNKS=()
for f in "${GENOME_DIR}"/chr25_part*_sub*[ab].fasta; do
    [[ -f "$f" ]] || continue
    name=$(basename "$f" .fasta)
    RESCUE_CHUNKS+=("$name")
done

if [[ ${#RESCUE_CHUNKS[@]} -eq 0 ]]; then
    echo "ERROR: No rescue sub-sub-chunk fastas found!"
    echo "Run resplit_chr25_rescue.py first."
    exit 1
fi

echo "Found ${#RESCUE_CHUNKS[@]} rescue chunks"

# ── Submit jobs ───────────────────────────────────────────────────
SUBMITTED=0
SKIPPED=0

for rescue in "${RESCUE_CHUNKS[@]}"; do
    # Parse: chr25_part01_sub06a → chr=chr25, parent=chr25_part01_sub06
    chr="${rescue%%_part*}"
    parent="${rescue%[ab]}"

    # Build job name: chr25_part01_sub06a → H25p01s06a
    part_num=$(echo "${parent}" | sed 's/chr25_part//' | cut -d'_' -f1)  # 01
    sub_num=$(echo "${parent}" | sed 's/chr25_part//' | cut -d'_' -f2 | sed 's/sub//')  # 06
    suffix="${rescue: -1}"  # a or b
    job_name="H25p${part_num}s${sub_num}${suffix}"

    # Output directory
    out_dir="${SCRIPT_DIR}/HiCAT_genome/${rescue}"

    # Skip if already completed
    if [[ -d "${out_dir}/out" ]] && ls "${out_dir}"/out_final_hor12.xls >/dev/null 2>&1; then
        echo "[skip] ${rescue} — already complete (layer 12)"
        ((SKIPPED++)) || true
        continue
    fi

    # Also skip if running/pending
    if squeue -u jhc103 -h -o "%j" 2>/dev/null | grep -q "^${job_name}$"; then
        echo "[skip] ${rescue} — already in queue as ${job_name}"
        ((SKIPPED++)) || true
        continue
    fi

    # Submit
    cmd="sbatch \
        --job-name=${job_name} \
        --mem=${MEM} \
        -c ${CPUS} \
        ${EXCLUDE_ARG} \
        ${SCRIPT_DIR}/run_HiCAT_rescue.sh ${rescue}"

    echo "[submit] ${job_name} ← ${rescue}  (${MEM}, ${CPUS}c)"

    if [[ "${DRY_RUN}" == "1" ]]; then
        echo "  [DRY RUN] ${cmd}"
        ((SUBMITTED++)) || true
    else
        result=$(eval "${cmd}")
        job_id=$(echo "${result}" | grep -oP '(?<=Submitted batch job )\d+')
        echo "  Job ID: ${job_id}"
        ((SUBMITTED++)) || true
    fi
done

echo ""
echo "=== Summary ==="
echo "Submitted: ${SUBMITTED}"
echo "Skipped:   ${SKIPPED}"
echo "Total:     ${#RESCUE_CHUNKS[@]}"
echo ""
echo "Monitor:   squeue -u jhc103 | grep H25p"
