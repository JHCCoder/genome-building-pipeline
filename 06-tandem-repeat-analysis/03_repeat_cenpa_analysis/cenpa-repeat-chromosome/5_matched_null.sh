#!/bin/bash
#SBATCH -J 0807_matched_null
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH --mem=16G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## Generate matched-shuffle null placements for each repeat array.
## For each array: sample N_SHUF null placements on the same chromosome with
## the exact same length, whose window covariates (mappability, GC, repeat
## density) fall within tolerance of the array's own window.
##
## Output: results/matched_null_placements.csv
##   chrom array_family array_start array_end start end mapp GC repdensity ok
## (ok=1 if a matched placement was found; null placements only when matched)

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh
init_dirs

THREADS=$SLURM_CPUS_PER_TASK
log "=== Step 5: matched shuffle null ==="

MAP="${MAPPABILITY_DIR}/mappability_win100kb.bed"
GC="${GC_DIR}/gc_win100kb.bed"
REPD="${REPDENS_DIR}/repdensity_win100kb.bed"
WINDOWS="${COVW_DIR}/windows_100kb.bed"
DOMAINS="${DOMAIN_DIR}/cenpa_domains_weighted.bed"
check_file "$MAP" "mappability"; check_file "$GC" "GC"; check_file "$REPD" "repeat density"
check_file "$DOMAINS" "domains"

python scripts/5a_matched_null.py \
    --arrays-bin6 "$BIN6_349" \
    --arrays-bin4 "$BIN4_195" \
    --arrays-bin8 "$BIN8_389" \
    --mapp "$MAP" --gc "$GC" --repd "$REPD" \
    --chrom-sizes "$CHROM_SIZES" \
    --n-shuf "$N_SHUFFLE" --seed "$SHUFFLE_SEED" \
    --out "${RESULTS_DIR}/matched_null_placements.csv"

log "=== Step 5 done ==="
