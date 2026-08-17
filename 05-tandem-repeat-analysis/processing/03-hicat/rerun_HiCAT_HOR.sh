#!/bin/bash
#SBATCH -J HiCAT_HOR_rerun
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4                    # low threads to avoid OOM in ED matrix calc
#SBATCH -t 4-00:00:00           # 4 days — ED ~7h (4 thr), pre-clustering N² ~3h, layers ~20h
#SBATCH --mem=300G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

# ============================================================
# Re-run ONLY the HiCAT_HOR.py step on existing decomposition.
#
# Use when: StringDecomposer completed but HiCAT_HOR crashed
#   silently (out_block.sequences exists but no out_all_layer*.xls).
#
# Usage: sbatch --job-name=HiCAT_rr_p03s01 rerun_HiCAT_HOR.sh chr4_part03_sub01
# ============================================================

source /tscc/nfs/home/jhc103/.bashrc

PART=${1:?Usage: sbatch rerun_HiCAT_HOR.sh <part_dir_name> [threads]}
HICAT_THREADS=${2:-4}  # default 4 threads; use 2 for dense regions (chr25)

BASE_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
SCRIPT_DIR="${BASE_DIR}/code/command-line-script/genome-annotation/HiCAT"
GENOME_DIR="${SCRIPT_DIR}/HiCAT_genome"
OUTPUT_DIR="${GENOME_DIR}/${PART}"

HICAT_HOR="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin/HiCAT_HOR.py"
DECOMP="${OUTPUT_DIR}/final_decomposition.tsv"
BASE_SEQ="${OUTPUT_DIR}/input_fasta.1.fa"
LOG="${OUTPUT_DIR}/log.txt"
BLOCK_SEQ="${OUTPUT_DIR}/out_block.sequences"

# ── Validate prerequisites ───────────────────────────────────────────
if [[ ! -f "${DECOMP}" ]]; then
    echo "ERROR: Decomposition not found: ${DECOMP}"
    echo "Run full HiCAT first (StringDecomposer step is missing)."
    exit 1
fi

if [[ ! -f "${BASE_SEQ}" ]]; then
    echo "ERROR: Base sequence not found: ${BASE_SEQ}"
    exit 1
fi

DECOMP_LINES=$(wc -l < "${DECOMP}")
if [[ "${DECOMP_LINES}" -le 1 ]]; then
    echo "ERROR: Decomposition is empty or header-only (${DECOMP_LINES} lines)"
    exit 1
fi

echo "[$(date)] ========================================"
echo "  Re-running HiCAT_HOR on: ${PART}"
echo "  Decomposition: ${DECOMP} (${DECOMP_LINES} lines)"
echo "  Base sequence: ${BASE_SEQ}"
echo "  Output dir:    ${OUTPUT_DIR}"
echo "  Threads:       ${HICAT_THREADS}"
echo "[$(date)] ========================================"

# ── Clean up partial output from previous failed run ─────────────────
# Remove files that HiCAT_HOR will regenerate, so we don't mistake
# stale partial output for a successful run.
rm -f "${OUTPUT_DIR}"/out_all_layer*.xls
rm -f "${OUTPUT_DIR}"/out_top_layer*.xls
rm -f "${OUTPUT_DIR}"/out_final_hor*.xls
rm -f "${OUTPUT_DIR}"/out_cluster_*.xls
rm -f "${OUTPUT_DIR}"/out_monomer_seq_*.xls
rm -f "${OUTPUT_DIR}"/out_pre_merge.xls
rm -f "${OUTPUT_DIR}"/out_statistics.xls
rm -f "${OUTPUT_DIR}"/out_block.sequences
rm -f "${OUTPUT_DIR}"/out_merge_edit_distance_matrix.xls
rm -f "${LOG}"
rm -rf "${OUTPUT_DIR}/out"

echo "Cleaned up partial files from previous run."

# ── Run HiCAT_HOR ────────────────────────────────────────────────────
export PATH="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin:$PATH"

python "${HICAT_HOR}" \
    -d "${DECOMP}" \
    -b "${BASE_SEQ}" \
    -o "${OUTPUT_DIR}" \
    -s 0.94 \
    -st 0.005 \
    -m 40 \
    -sp 5 \
    -sn 10 \
    -t ${HICAT_THREADS}

EXIT_CODE=$?
echo "[$(date)] HiCAT_HOR finished with exit code ${EXIT_CODE}"

# ── Validate output ─────────────────────────────────────────────────
LAYER_COUNT=$(ls "${OUTPUT_DIR}"/out_all_layer*.xls 2>/dev/null | wc -l)

echo ""
echo "Layer files produced: ${LAYER_COUNT}"
ls -lh "${OUTPUT_DIR}"/out_all_layer*.xls 2>/dev/null
echo ""

if [[ -d "${OUTPUT_DIR}/out" ]]; then
    echo "Final out/ dir:"
    ls -lh "${OUTPUT_DIR}/out/"
else
    echo "WARNING: out/ directory not created (final HOR selection step)"
fi

if [[ "${LAYER_COUNT}" -eq 0 ]]; then
    echo ""
    echo "ERROR: No layer files produced. HiCAT_HOR failed silently."
    echo "Last log entry:"
    cat "${LOG}" 2>/dev/null || echo "  (no log file)"
    exit 1
fi

echo ""
echo "SUCCESS: ${PART} HTRM clustering complete."
exit 0
