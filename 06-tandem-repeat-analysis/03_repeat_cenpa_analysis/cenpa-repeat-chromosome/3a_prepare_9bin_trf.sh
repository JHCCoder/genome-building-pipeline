#!/bin/bash
#SBATCH -J 0807_9bin_prep
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## Build merged (per-bin) TRF BEDs from the period-binned TRF file, and a
## merged "all-repeats" BED for the repeat-density track.
## TRF bin BED format: chr start end interval_id bin_id period_size copies match_pct

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh
init_dirs

REPDENS_DIR="${REPDENS_DIR}"   # from config
mkdir -p "${REPDENS_DIR}/per_bin"

# per-bin merged beds
for b in 1 2 3 4 5 6 7 8 9; do
    out="${REPDENS_DIR}/per_bin/bin${b}_merged.bed"
    if [[ -f "$out" ]] && [[ -s "$out" ]]; then continue; fi
    awk -v b="$b" '$5==b' "$TRF_BED" | cut -f1-3 \
        | sort -k1,1V -k2,2n \
        | bedtools merge -d 0 -i - > "$out"
    log "bin$b: $(wc -l < "$out") merged loci"
done

# all-repeats (any bin) union
ALL_OUT="${REPDENS_DIR}/all_repeats_merged.bed"
if [[ ! -s "$ALL_OUT" ]]; then
    cut -f1-3 "$TRF_BED" | sort -k1,1V -k2,2n \
        | bedtools merge -d 0 -i - > "$ALL_OUT"
    log "all-repeats merged: $(wc -l < "$ALL_OUT") loci"
fi

log "=== 3a done ==="
