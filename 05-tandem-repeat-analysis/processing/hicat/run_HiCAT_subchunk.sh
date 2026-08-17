#!/bin/bash
#SBATCH -J HiCAT_sub           # override with --job-name=HiCAT_c4p01s01
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4                    # threads — keep low: HiCAT_HOR ED matrix is N² memory; 32 threads = 32 worker copies
#SBATCH -t 4-00:00:00           # 4 days — fewer threads = slower ED calc, but avoids OOM
#SBATCH --mem=300G              # peak ~200 GB for N≈90k (3 matrix copies); 300G gives headroom
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type FAIL
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

# ============================================================
# HiCAT sub-chunk runner (smaller pieces of failed chunks)
#
# Usage: sbatch --job-name=HiCAT_c4p01s01 run_HiCAT_subchunk.sh chr4_part01_sub01
#
# Naming convention:
#   chr4_part01_sub01  →  chromosome chr4, original part 01, sub-chunk 01
#   Fasta header fixed: chr4_part01_sub01 → chr4p01s01
#
# Provenance: mergable back to chr4_part01 by concatenating
#   all chr4_part01_sub*/final_decomposition.tsv
# ============================================================

source /tscc/nfs/home/jhc103/.bashrc

PART=${1:?Usage: sbatch run_HiCAT_subchunk.sh <chr_part_sub, e.g. chr4_part01_sub01>}

BASE_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
SCRIPT_DIR="${BASE_DIR}/code/command-line-script/genome-annotation/HiCAT"
TEMPLATE_DIR="${BASE_DIR}/output/outputs-from-centraAnno/hifiasm-0414"
GENOME_DIR="${SCRIPT_DIR}/HiCAT_genome"

# ── Parse chromosome and original part from sub-chunk name ──────────
# chr4_part01_sub01 → CHR=chr4, ORIG_PART=chr4_part01, SUB=01
CHR="${PART%%_part*}"
ORIG_PART="${PART%_sub*}"          # chr4_part01_sub01 → chr4_part01
SUB_NUM="${PART##*_sub}"           # chr4_part01_sub01 → 01

CHUNK_FASTA="${GENOME_DIR}/${CHR}/${PART}.fasta"
OUTPUT_DIR="${GENOME_DIR}/${PART}"

# ── Select centromere-specific template ─────────────────────────────
case "${CHR}" in
    chr4)
        TEMPLATE="${TEMPLATE_DIR}/chr4_monomerTemplates_line1743-1804.fa"
        ;;
    chr25)
        TEMPLATE="${TEMPLATE_DIR}/chr25_monomerTemplates_centromere.fa"
        ;;
    *)
        echo "ERROR: No centromere-specific template defined for ${CHR}"
        exit 1
        ;;
esac

# ── Validate inputs ─────────────────────────────────────────────────
if [[ ! -f "${TEMPLATE}" ]]; then
    echo "ERROR: Template not found: ${TEMPLATE}"
    exit 1
fi

if [[ ! -f "${CHUNK_FASTA}" ]]; then
    echo "ERROR: Sub-chunk fasta not found: ${CHUNK_FASTA}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# ── Fix fasta header ────────────────────────────────────────────────
# HiCAT_HOR.py splits monomer names on "_" and expects integer coords.
# chr4_part01_sub01 → chr4p01s01 (no underscores, preserves provenance)
ORIG_PART_NUM="${ORIG_PART##*_part}"   # chr4_part01 → 01
FIXED_NAME="${CHR}p${ORIG_PART_NUM}s${SUB_NUM}"

FIXED_FASTA="${OUTPUT_DIR}/input_fixed.fa"
sed "s/>${PART}/>${FIXED_NAME}/" "${CHUNK_FASTA}" > "${FIXED_FASTA}"

echo "[$(date)] ========================================"
echo "  Sub-chunk:   ${PART}"
echo "  Parent:      ${ORIG_PART}"
echo "  Chromosome:  ${CHR}"
echo "  Fixed name:  ${FIXED_NAME}"
echo "  Input:       ${FIXED_FASTA} ($(wc -c < ${FIXED_FASTA}) bp)"
echo "  Template:    ${TEMPLATE} ($(grep -c '^>' ${TEMPLATE}) monomers)"
echo "  Output:      ${OUTPUT_DIR}"
echo "  Threads:     ${SLURM_CPUS_PER_TASK:-16}"
echo "  Memory:      ${SLURM_MEM_PER_NODE:-128G}"
echo "[$(date)] ========================================"

# ── Run HiCAT ───────────────────────────────────────────────────────
export PATH="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin:$PATH"

# Use 4 threads — HiCAT_HOR ED matrix is O(N²) memory and joblib spawns
# N workers, each holding a copy of the data. 32 workers = OOM; 4 = safe.
HICAT_THREADS=4

/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin/hicat \
    --output_dir "${OUTPUT_DIR}" \
    -i "${FIXED_FASTA}" \
    -t "${TEMPLATE}" \
    -th "${HICAT_THREADS}"

EXIT_CODE=$?
echo "[$(date)] HiCAT finished for ${PART} with exit code ${EXIT_CODE}"

# ── Validate output ─────────────────────────────────────────────────
# HiCAT_HOR.py can crash silently (hicat wrapper exits 0). Check that
# the final output was actually produced.
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
    echo "PROVENANCE: ${ORIG_PART} ← ${PART}"

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
