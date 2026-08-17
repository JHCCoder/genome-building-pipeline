#!/bin/bash
#SBATCH -J HiCAT_chr           # override with --job-name=HiCAT_chr1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32                   # threads for HiCAT
#SBATCH -t 4-00:00:00           # default 4 days; override per-chr below
#SBATCH --mem=250G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

# ============================================================
# HiCAT per-chromosome runner
# Usage: sbatch [--job-name=HiCAT_chr1] [--mem=300G] [-t 7-00:00:00] run_HiCAT_chr.sh <chr>
# ============================================================

source /tscc/nfs/home/jhc103/.bashrc

CHR=${1:?Usage: sbatch run_HiCAT_chr.sh <chromosome, e.g. chr1>}

BASE_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
ASSEMBLY="${BASE_DIR}/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta"
TEMPLATE_DIR="${BASE_DIR}/output/outputs-from-centraAnno/hifiasm-0414"
OUTPUT_DIR="${BASE_DIR}/code/command-line-script/genome-annotation/HiCAT/HiCAT_genome/${CHR}"

mkdir -p "${OUTPUT_DIR}"

# --------------------------------------------------
# Extract chromosome FASTA from full assembly
# --------------------------------------------------
CHR_FASTA="${OUTPUT_DIR}/${CHR}.fasta"
echo "[$(date)] Extracting ${CHR} from assembly ..."

/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-samtools/bin/samtools faidx "${ASSEMBLY}" "${CHR}" | \
	    awk '/^>/ {print; next} {print toupper($0)}' > "${CHR_FASTA}"

if [[ ! -s "${CHR_FASTA}" ]]; then
    echo "ERROR: Failed to extract ${CHR} from assembly"
    exit 1
fi

echo "[$(date)] ${CHR} extracted: $(wc -c < ${CHR_FASTA}) bp"

# --------------------------------------------------
# Run HiCAT
# --------------------------------------------------
TEMPLATE="${TEMPLATE_DIR}/${CHR}_monomerTemplates.fa"

if [[ ! -f "${TEMPLATE}" ]]; then
    echo "ERROR: Monomer template not found: ${TEMPLATE}"
    exit 1
fi

echo "[$(date)] Starting HiCAT for ${CHR} ..."
echo "  Input:    ${CHR_FASTA}"
echo "  Template: ${TEMPLATE}"
echo "  Output:   ${OUTPUT_DIR}"
echo "  Threads:  ${SLURM_CPUS_PER_TASK:-32}"

export PATH="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin:$PATH"
/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin/hicat \
    --output_dir "${OUTPUT_DIR}" \
    -i "${CHR_FASTA}" \
    -t "${TEMPLATE}" \
    -th "${SLURM_CPUS_PER_TASK:-32}"

EXIT_CODE=$?
echo "[$(date)] HiCAT finished for ${CHR} with exit code ${EXIT_CODE}"

exit ${EXIT_CODE}
