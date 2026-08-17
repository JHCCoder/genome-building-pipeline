#!/bin/bash
# ============================================================
# Submit all sub-chunk HiCAT jobs with provenance tracking.
#
# Sub-chunks are created by resplit_failed_chunks.py.
# Resource baseline from job 7550536:
#   Full 10Mb chr4 centromere, ~86k blocks, 300G, 8 cores
#   Completed 18h, 159 GB actual usage, 53% memory efficiency.
#
# All sub-chunks target ~90k blocks → 300G, 32 cores, 2d.
#
# Manifest: HiCAT_genome/subchunk_manifest.tsv
#
# Usage: bash submit_subchunks.sh [--dry-run]
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GENOME_DIR="${SCRIPT_DIR}/HiCAT_genome"
MANIFEST="${GENOME_DIR}/subchunk_manifest.tsv"

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

# ── All sub-chunks: same resources as reference job 7550536 ──────────
# Targeting ~90k blocks → 300G, 32 cores
MEM="300G"
CPU=32

# ── Discover all sub-chunk fastas ───────────────────────────────────
echo "=== Discovering sub-chunk fastas ==="
SUBCHUNKS=()
for f in "${GENOME_DIR}"/chr4/chr4_part*_sub*.fasta \
         "${GENOME_DIR}"/chr25/chr25_part*_sub*.fasta; do
    [[ -f "$f" ]] || continue
    name=$(basename "$f" .fasta)         # chr4_part01_sub01
    SUBCHUNKS+=("$name")
done

if [[ ${#SUBCHUNKS[@]} -eq 0 ]]; then
    echo "ERROR: No sub-chunk fastas found!"
    echo "Run resplit_failed_chunks.py first."
    exit 1
fi

echo "Found ${#SUBCHUNKS[@]} sub-chunks"

# ── Write manifest header ───────────────────────────────────────────
if [[ "${DRY_RUN}" != "1" ]]; then
    echo -e "sub_chunk\tchromosome\tparent_chunk\tsub_index\tjob_id\tmem\tcpus\tfasta_path" > "${MANIFEST}"
fi

# ── Submit jobs ─────────────────────────────────────────────────────
SUBMITTED=0
SKIPPED=0

for sub in "${SUBCHUNKS[@]}"; do
    # Parse: chr4_part01_sub01 → chr=chr4, parent=chr4_part01, sub_idx=01
    chr="${sub%%_part*}"
    parent="${sub%_sub*}"
    sub_idx="${sub##*_sub}"

    # Job name: HiCAT_c4p01s01
    chr_short="${chr/chr/c}"           # chr4 → c4
    parent_num="${parent##*_part}"     # chr4_part01 → 01
    job_name="HiCAT_${chr_short}p${parent_num}s${sub_idx}"

    # Output directory
    out_dir="${GENOME_DIR}/${sub}"

    # Skip if already completed (has out_ files)
    if ls "${out_dir}"/out_* >/dev/null 2>&1; then
        echo "[skip] ${sub} — already has output files"
        ((SKIPPED++)) || true
        continue
    fi

    # Submit
    cmd="sbatch \
        --job-name=${job_name} \
        --mem=${MEM} \
        -c ${CPU} \
        ${SCRIPT_DIR}/run_HiCAT_subchunk.sh ${sub}"

    echo "[submit] ${job_name} ← ${sub}  (${MEM}, ${CPU}c)"

    if [[ "${DRY_RUN}" == "1" ]]; then
        echo "  [DRY RUN] ${cmd}"
        job_id="DRY_RUN"
    else
        result=$(eval "${cmd}")
        job_id=$(echo "${result}" | grep -oP '(?<=Submitted batch job )\d+')
        echo "  Job ID: ${job_id}"
    fi

    # Record in manifest
    if [[ "${DRY_RUN}" != "1" ]]; then
        echo -e "${sub}\t${chr}\t${parent}\t${sub_idx}\t${job_id}\t${MEM}\t${CPU}\t${GENOME_DIR}/${chr}/${sub}.fasta" >> "${MANIFEST}"
    fi

    ((SUBMITTED++)) || true
done

echo ""
echo "=== Summary ==="
echo "Submitted: ${SUBMITTED}"
echo "Skipped:   ${SKIPPED}"
echo "Total:     ${#SUBCHUNKS[@]}"
echo ""
echo "Manifest:  ${MANIFEST}"
echo ""
echo "Monitor:   squeue -u jhc103 | grep HiCAT_c"
echo ""
echo "Provenance columns in manifest:"
echo "  sub_chunk    → parent_chunk (merge key)"
echo "  e.g. chr4_part01_sub01, chr4_part01_sub02 → chr4_part01"
