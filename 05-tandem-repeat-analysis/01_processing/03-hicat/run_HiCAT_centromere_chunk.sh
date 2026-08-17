#!/bin/bash
#SBATCH -J HiCAT_cen           # override with --job-name
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32                   # default; override per-chunk below
#SBATCH -t 2-00:00:00           # default 2 days; override per-chunk below
#SBATCH --mem=250G              # default; override per-chunk below
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

# ============================================================
# HiCAT centromere-specific chunk runner
# Usage: sbatch run_HiCAT_centromere_chunk.sh <chromosome> <part_number>
#
# Uses centromere-specific monomer templates (identified from
# the largest HOR arrays per chromosome) to keep monomer counts
# manageable for HiCAT_HOR.py's O(n^2) algorithm.
#
# Fixes chunk fasta header naming: chr4_part01 → chr4p01
# to avoid underscore-parsing crash in HiCAT_HOR.py.
# ============================================================

source /tscc/nfs/home/jhc103/.bashrc

CHR=${1:?Usage: sbatch run_HiCAT_centromere_chunk.sh <chr, e.g. chr4> <part, e.g. 01>}
PART_NUM=${2:?Usage: sbatch run_HiCAT_centromere_chunk.sh <chr, e.g. chr4> <part, e.g. 01>}

BASE_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
SCRIPT_DIR="${BASE_DIR}/code/command-line-script/genome-annotation/HiCAT"
TEMPLATE_DIR="${BASE_DIR}/output/outputs-from-centraAnno/hifiasm-0414"
GENOME_DIR="${SCRIPT_DIR}/HiCAT_genome"

PART="${CHR}_part${PART_NUM}"
CHUNK_FASTA="${GENOME_DIR}/${CHR}/${PART}.fasta"
OUTPUT_DIR="${GENOME_DIR}/${PART}"

# --------------------------------------------------
# Select centromere-specific template per chromosome
# --------------------------------------------------
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

if [[ ! -f "${TEMPLATE}" ]]; then
    echo "ERROR: Template not found: ${TEMPLATE}"
    exit 1
fi

if [[ ! -f "${CHUNK_FASTA}" ]]; then
    echo "ERROR: Chunk fasta not found: ${CHUNK_FASTA}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# --------------------------------------------------
# Fix fasta header: chr4_part01 → chr4p01
# HiCAT_HOR.py splits monomer names on "_" and expects
# integer coordinates, so sequence names must not contain
# underscores (chr4_part01 → chr4p01).
# --------------------------------------------------
FIXED_NAME="${CHR}p${PART_NUM}"
FIXED_FASTA="${OUTPUT_DIR}/input_fixed.fa"

sed "s/>${PART}/>${FIXED_NAME}/" "${CHUNK_FASTA}" > "${FIXED_FASTA}"

echo "[$(date)] ========================================"
echo "  Chromosome:  ${CHR}"
echo "  Part:        ${PART}"
echo "  Fixed name:  ${FIXED_NAME}"
echo "  Input:       ${FIXED_FASTA} ($(wc -c < ${FIXED_FASTA}) bp)"
echo "  Template:    ${TEMPLATE} ($(grep -c '^>' ${TEMPLATE}) monomers)"
echo "  Output:      ${OUTPUT_DIR}"
echo "  Threads:     ${SLURM_CPUS_PER_TASK:-32}"
echo "  Memory:      ${SLURM_MEM_PER_NODE:-250G}"
echo "[$(date)] ========================================"

# --------------------------------------------------
# Run HiCAT
# --------------------------------------------------
export PATH="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin:$PATH"

/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin/hicat \
    --output_dir "${OUTPUT_DIR}" \
    -i "${FIXED_FASTA}" \
    -t "${TEMPLATE}" \
    -th "${SLURM_CPUS_PER_TASK:-32}"

EXIT_CODE=$?
echo "[$(date)] HiCAT finished for ${PART} with exit code ${EXIT_CODE}"

# --------------------------------------------------
# Report output files
# --------------------------------------------------
if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo ""
    echo "Output files:"
    ls -lh "${OUTPUT_DIR}"/out_* 2>/dev/null || echo "  (no out_* files yet — check log)"
    echo ""
    echo "Decomposition:"
    wc -l "${OUTPUT_DIR}"/final_decomposition.tsv 2>/dev/null || echo "  (not found)"
fi

exit ${EXIT_CODE}
