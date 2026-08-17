#!/bin/bash
#SBATCH -J 0807_nh_weighted
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 8:00:00
#SBATCH --mem=64G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## Usage: sbatch --array=1-3 1_nh_weighted.sh   (150,151,152)

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh
init_dirs

SAMPLES=(XG_150 XG_151 XG_152)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE=${SAMPLES[$IDX]}

case "$SAMPLE" in
  XG_150) BAM="$BAM_150_K100"; FRAG="$FRAG_150" ;;
  XG_151) BAM="$BAM_151_K100"; FRAG="$FRAG_151" ;;
  XG_152) BAM="$BAM_152_K100"; FRAG="$FRAG_152" ;;
esac

THREADS=$SLURM_CPUS_PER_TASK
log "=== Step 1: $SAMPLE 1/NH weighted coverage ==="
check_file "$BAM" "k100 BAM"
check_file "$FRAG" "k1 fragment BED"

# --- k=100 weighted window coverage (1/NH) ---
python scripts/1a_compute_nh_weighted.py "$SAMPLE" "$BAM" \
    "$CHROM_SIZES" "$COVW_DIR" "$THREADS" 2>&1 | tee "${COVW_DIR}/${SAMPLE}_step1.log"

# --- k=1 unique-mapping window coverage (sensitivity leg) ---
WINDOWS="${COVW_DIR}/windows_100kb.bed"
bedtools coverage -a "$WINDOWS" -b "$FRAG" -counts \
    > "${COVW_DIR}/${SAMPLE}_win100kb_k1.txt"

log "=== Step 1 done: $SAMPLE ==="
