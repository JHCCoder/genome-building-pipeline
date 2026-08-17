#!/bin/bash
#SBATCH -J 0807_define_probes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH --mem=4G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --no-requeue

## Select family-specific 31-mer probes from the 6c count files (fast).

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh

python scripts/6c_define_probes.py "$PROBE_DIR" 2>&1 | tee "${PROBE_DIR}/probe_definition.log"
log "=== probe definition done ==="
