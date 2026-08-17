#!/bin/bash
#SBATCH --job-name=081225_telomereTRF_thorough
#SBATCH --output=/tscc/nfs/home/jhc103/cluster-logs/081225_telomereTRF_thorough.%A.%a.%N.out
#SBATCH --error=/tscc/nfs/home/jhc103/cluster-logs/081225_telomereTRF_thorough.%A.%a.%N.err
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBATCH --account=csd788
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=1-30

# Thorough genome-wide TRF scan for telomeric tandem repeats.
# More sensitive than the original 2 7 7 run: lower mismatch penalty (5 vs 7)
# catches degenerate telomere arrays (imperfect TTAGGG / CCCTAA repeats).
#
# One array task per chromosome (chr1-28, chrX, chrY).

set -euo pipefail

TRF="/tscc/projects/ps-renlab2/jhc103/toolshed/repeat-annotation/trf409.linux64"
CHROM_DIR="/tscc/lustre/ddn/scratch/jhc103/telomere-thorough-scan/chroms"
OUT_DIR="/tscc/lustre/ddn/scratch/jhc103/telomere-thorough-scan/trf_out"

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
        chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chrX chrY)

i=$(( SLURM_ARRAY_TASK_ID - 1 ))
chrom="${CHROMS[$i]}"
input="${CHROM_DIR}/${chrom}.fa"
cd "$OUT_DIR" || exit 1

echo "[$(date '+%F %T')] Task ${SLURM_ARRAY_TASK_ID}/30: TRF on ${chrom}"

"$TRF" "$input" 2 5 7 80 10 50 500 -m -d -l 6

echo "[$(date '+%F %T')] Done ${chrom}"

# Data file (chrX.fa.2.5.7.80.10.50.500.dat) now exists in OUT_DIR
ls -la "${chrom}.fa.2.5.7.80.10.50.500.dat"
