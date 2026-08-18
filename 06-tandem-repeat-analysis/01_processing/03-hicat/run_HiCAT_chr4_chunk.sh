#!/bin/bash
#SBATCH -J HiCAT_chunk       # override with --job-name=HiCAT_c4_p01, etc.
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4                    # threads — keep low to avoid OOM in HiCAT_HOR ED matrix calc
#SBATCH -t 4-00:00:00           # 4 days — fewer threads = slower but won't OOM
#SBATCH --mem=300G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

# ============================================================
# HiCAT chromosome chunk runner (general-purpose)
# Usage: sbatch run_HiCAT_chunk.sh <chr_part, e.g. chr4_part01, chr25_part03>
#
# The part name encodes both chromosome and chunk: chrN_partMM
# Chunk fastas must be pre-split via split_chromosome.py into
#   HiCAT_genome/<chr>/<chr>_partMM.fasta
# ============================================================

source /tscc/nfs/home/jhc103/.bashrc

PART=${1:?Usage: sbatch run_HiCAT_chunk.sh <chr_part, e.g. chr25_part01>}

BASE_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
SCRIPT_DIR="${BASE_DIR}/code/command-line-script/genome-annotation/HiCAT"
TEMPLATE_DIR="${BASE_DIR}/output/outputs-from-centraAnno/hifiasm-0414"
OUTPUT_DIR="${SCRIPT_DIR}/HiCAT_genome/${PART}"

mkdir -p "${OUTPUT_DIR}"

# --------------------------------------------------
# Map part → chromosome for template + input lookup
# e.g. chr25_part03 → CHR=chr25
# --------------------------------------------------
CHR="${PART%%_part*}"

# Locate input fasta (pre-split chunk, in HiCAT_genome/<chr>/)
CHUNK_FASTA="${SCRIPT_DIR}/HiCAT_genome/${CHR}/${PART}.fasta"

if [[ ! -s "${CHUNK_FASTA}" ]]; then
    echo "ERROR: Chunk fasta not found: ${CHUNK_FASTA}"
    exit 1
fi

echo "[$(date)] HiCAT starting for ${PART}"
echo "  Input:    ${CHUNK_FASTA}"
echo "  Size:     $(wc -c < ${CHUNK_FASTA}) bp"

# Template lookup (also uses CHR extracted above)
TEMPLATE="${TEMPLATE_DIR}/${CHR}_monomerTemplates.fa"

if [[ ! -f "${TEMPLATE}" ]]; then
    echo "ERROR: Monomer template not found: ${TEMPLATE}"
    exit 1
fi

# ── Fix fasta header ────────────────────────────────────────────────
# HiCAT_HOR.py splits monomer names on "_" and expects integer coords.
# Names like "chr4_part01" break the parsing → ValueError crash.
# Replace underscores: chr4_part01 → chr4p01
CHR_NUM="${PART##*_part}"
FIXED_NAME="${CHR}p${CHR_NUM}"
FIXED_FASTA="${OUTPUT_DIR}/input_fixed.fa"
sed "s/>${PART}/>${FIXED_NAME}/" "${CHUNK_FASTA}" > "${FIXED_FASTA}"

echo "  Template: ${TEMPLATE}"
echo "  Output:   ${OUTPUT_DIR}"
echo "  Fixed name: ${FIXED_NAME}"
echo "  Threads:  4"

# --------------------------------------------------
# Run HiCAT (4 threads — higher values cause OOM in ED matrix calc)
# --------------------------------------------------
export PATH="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin:$PATH"
/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin/hicat \
    --output_dir "${OUTPUT_DIR}" \
    -i "${FIXED_FASTA}" \
    -t "${TEMPLATE}" \
    -th 4

EXIT_CODE=$?
echo "[$(date)] HiCAT finished for ${PART} with exit code ${EXIT_CODE}"

# ── Validate output ─────────────────────────────────────────────────
if [[ ${EXIT_CODE} -eq 0 ]]; then
    LAYER_COUNT=$(ls "${OUTPUT_DIR}"/out_all_layer*.xls 2>/dev/null | wc -l)
    if [[ "${LAYER_COUNT}" -eq 0 ]]; then
        echo "ERROR: HiCAT exited 0 but produced no HTRM output (silent crash)"
        exit 1
    fi
fi

exit ${EXIT_CODE}
