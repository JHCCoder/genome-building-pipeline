#!/bin/bash
#SBATCH -J HiCAT_rescue        # override with --job-name
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2                    # 2 threads — conservative, less memory pressure
#SBATCH -t 2-00:00:00           # 2 days — 47k blocks should finish in <24h
#SBATCH --mem=300G              # 300G — ED matrix ~18 GB/copy × 3 copies ≈ 54 GB; plenty of headroom
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type FAIL
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

# ============================================================
# HiCAT chr25 rescue runner — 2-way split sub-sub-chunks
#
# Usage: sbatch --job-name=H25p01s06a run_HiCAT_rescue.sh chr25_part01_sub06a
#
# Naming convention:
#   chr25_part01_sub06a  →  chr25 part 01, sub-chunk 06, rescue half a
#   Fasta header fixed: chr25_part01_sub06a → chr25p01s06a
#
# Provenance: mergable back to chr25_part01_sub06 by
#   concatenating a and b final_decomposition.tsv
#
# Conservative settings (vs original 500G/4 threads):
#   - 300G memory (5× actual need for 47k blocks)
#   - 2 threads (fewer ED matrix copies → less memory)
#   - GPU nodes excluded (cluster instability on GPU nodes)
# ============================================================

source /tscc/nfs/home/jhc103/.bashrc

PART=${1:?Usage: sbatch run_HiCAT_rescue.sh <chr25_partNN_subMM[a|b]>}

BASE_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
SCRIPT_DIR="${BASE_DIR}/code/command-line-script/genome-annotation/HiCAT"
TEMPLATE_DIR="${BASE_DIR}/output/outputs-from-centraAnno/hifiasm-0414"
GENOME_DIR="${SCRIPT_DIR}/HiCAT_genome"

# ── Parse chromosome, parent sub-chunk, and rescue suffix ──────────
# chr25_part01_sub06a → CHR=chr25, PARENT=chr25_part01_sub06, SUFFIX=a
CHR="${PART%%_part*}"
PARENT="${PART%[ab]}"                          # chr25_part01_sub06a → chr25_part01_sub06
SUFFIX="${PART: -1}"                           # last char: a or b

CHUNK_FASTA="${GENOME_DIR}/${CHR}/${PART}.fasta"
OUTPUT_DIR="${GENOME_DIR}/${PART}"

# ── Centromere-specific template ──────────────────────────────────
TEMPLATE="${TEMPLATE_DIR}/chr25_monomerTemplates_centromere.fa"

if [[ ! -f "${TEMPLATE}" ]]; then
    echo "ERROR: Template not found: ${TEMPLATE}"
    exit 1
fi

if [[ ! -f "${CHUNK_FASTA}" ]]; then
    echo "ERROR: Sub-sub-chunk fasta not found: ${CHUNK_FASTA}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# ── Fix fasta header ──────────────────────────────────────────────
# HiCAT_HOR.py splits monomer names on "_" and expects integer coords.
# chr25_part01_sub06a → chr25p01s06a (no underscores, preserves provenance)
PART_NUM="${PARENT#chr25_part}"        # 01_sub06
PART_NUM="${PART_NUM//_/}"             # 01sub06 → but we need p01s06 format
# Reconstruct: chr25_part01_sub06 → p01s06 + suffix
P=$(echo "${PARENT}" | sed 's/chr25_part//' | sed 's/_sub//')  # 0106
FIXED_NAME="chr25p${P:0:2}s${P:2:2}${SUFFIX}"                  # chr25p01s06a

FIXED_FASTA="${OUTPUT_DIR}/input_fixed.fa"
sed "s/>${PART}/>${FIXED_NAME}/" "${CHUNK_FASTA}" > "${FIXED_FASTA}"

echo "[$(date)] ========================================"
echo "  Rescue:      ${PART}"
echo "  Parent:      ${PARENT}"
echo "  Chromosome:  ${CHR}"
echo "  Suffix:      ${SUFFIX}"
echo "  Fixed name:  ${FIXED_NAME}"
echo "  Input:       ${FIXED_FASTA} ($(wc -c < ${FIXED_FASTA}) bp)"
echo "  Template:    ${TEMPLATE} ($(grep -c '^>' ${TEMPLATE}) monomers)"
echo "  Output:      ${OUTPUT_DIR}"
echo "  Threads:     2"
echo "  Memory:      300G"
echo "[$(date)] ========================================"

# ── Run HiCAT with conservative thread count ──────────────────────
export PATH="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin:$PATH"

# 2 threads — ED matrix is O(N²) memory per worker.
# At 47k blocks: ~18 GB per matrix copy × 3 copies ≈ 54 GB. 300G is comfortable.
HICAT_THREADS=2

/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin/hicat \
    --output_dir "${OUTPUT_DIR}" \
    -i "${FIXED_FASTA}" \
    -t "${TEMPLATE}" \
    -th "${HICAT_THREADS}"

EXIT_CODE=$?
echo "[$(date)] HiCAT finished for ${PART} with exit code ${EXIT_CODE}"

# ── Validate output ───────────────────────────────────────────────
if [[ ${EXIT_CODE} -eq 0 ]]; then
    LAYER_COUNT=$(ls "${OUTPUT_DIR}"/out_all_layer*.xls 2>/dev/null | wc -l)
    HAS_OUT_DIR=$(test -d "${OUTPUT_DIR}/out" && echo "yes" || echo "no")

    echo ""
    echo "Output files:"
    ls -lh "${OUTPUT_DIR}"/out_* 2>/dev/null || echo "  (no out_* files)"
    echo ""
    echo "Decomposition:"
    wc -l "${OUTPUT_DIR}"/final_decomposition.tsv 2>/dev/null || echo "  (not found)"
    echo ""
    echo "Layer files: ${LAYER_COUNT}"
    echo "Final out/ dir: ${HAS_OUT_DIR}"
    echo "PROVENANCE: ${PARENT} ← ${PART}"

    # Fail if no HTRM output was produced (silent crash detection)
    if [[ "${LAYER_COUNT}" -eq 0 ]]; then
        echo ""
        echo "ERROR: HiCAT exited 0 but produced no out_all_layer*.xls files."
        echo "This indicates a silent crash in HiCAT_HOR.py (likely memory)."
        echo "Check ${OUTPUT_DIR}/log.txt for last progress message."
        exit 1
    fi
    if [[ "${HAS_OUT_DIR}" == "no" ]]; then
        echo ""
        echo "ERROR: HiCAT exited 0 but out/ directory was not created."
        echo "The HTRM clustering did not complete."
        exit 1
    fi
fi

exit ${EXIT_CODE}
