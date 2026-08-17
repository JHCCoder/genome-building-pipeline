#!/bin/bash
#SBATCH -J 0807_covariates
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## Per-window covariates for the matched-shuffle null:
##   mappability   = fraction of 100kb window bp covered by >=1 uniquely-mapped
##                   control (H3K27ac) fragment  [empirical short-read mappability]
##   GC content    = (G+C) / (A+C+G+T) over the window from the assembly
##   repeat density= fraction of window bp covered by ANY merged TRF repeat
##
## NOTE: fragment BEDs include unplaced seq* contigs; restrict to chr1-28,X,Y
## and sort so bedtools coverage runs in streaming (-sorted) mode.

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

cd /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome
source scripts/config.sh
init_dirs

THREADS=$SLURM_CPUS_PER_TASK
log "=== Step 3: per-window covariates ==="

WINDOWS="${COVW_DIR}/windows_100kb.bed"
check_file "$WINDOWS" "windows"

# ---- mappability: pooled control (XG_152 + XG_153) fragment coverage ----
MAPOUT="${MAPPABILITY_DIR}/mappability_win100kb.bed"
if [[ ! -s "$MAPOUT" ]]; then
    POOLED="${MAPPABILITY_DIR}/pooled_ctrl_fragments_chr.bed"
    if [[ ! -s "$POOLED" ]]; then
        cat "$FRAG_152" "$FRAG_153" \
            | awk '$1 ~ /^chr([0-9]+|X|Y)$/' \
            | sort -k1,1V -k2,2n -T /tscc/lustre/ddn/scratch/jhc103 \
            > "$POOLED"
        log "pooled chr-only control fragments: $(wc -l < "$POOLED")"
    fi
    bedtools coverage -a "$WINDOWS" -b "$POOLED" -sorted \
        > "${MAPPABILITY_DIR}/mappability_win100kb_raw.bed"
    # bedtools coverage cols: chr start end n_features overlap_bp window_bp frac
    # col7 = fraction of window bp with >=1 fragment = mappability proxy
    awk -v OFS='\t' '{print $1, $2, $3, ($7==0 ? 0 : $7)}' \
        "${MAPPABILITY_DIR}/mappability_win100kb_raw.bed" > "$MAPOUT"
    log "mappability: $(wc -l < "$MAPOUT") windows"
fi

# ---- GC content (streaming; -nuc reads only the requested windows) ----
GCOUT="${GC_DIR}/gc_win100kb.bed"
if [[ ! -s "$GCOUT" ]]; then
    bedtools nuc -fi "$ASSEMBLY_FASTA" -bed "$WINDOWS" \
        > "${GC_DIR}/gc_win100kb_raw.bed"
    # nuc cols: chr start end pct_at pct_gc num_A num_C num_G num_T num_N num_oth seq_len
    awk -v OFS='\t' 'NR==1{next} {print $1, $2, $3, ($5==-1 ? 0 : $5)}' \
        "${GC_DIR}/gc_win100kb_raw.bed" > "$GCOUT"
    log "GC: $(wc -l < "$GCOUT") windows"
fi

# ---- repeat density (all TRF merged, chr-only) ----
REPOUT="${REPDENS_DIR}/repdensity_win100kb.bed"
ALLREP="${REPDENS_DIR}/all_repeats_merged.bed"
if [[ ! -s "$REPOUT" ]]; then
    check_file "$ALLREP" "all-repeats merged (run 3a first)"
    bedtools coverage -a "$WINDOWS" -b "$ALLREP" -sorted \
        > "${REPDENS_DIR}/repdensity_raw.bed"
    awk -v OFS='\t' '{print $1, $2, $3, ($7==0 ? 0 : $7)}' \
        "${REPDENS_DIR}/repdensity_raw.bed" > "$REPOUT"
    log "repeat density: $(wc -l < "$REPOUT") windows"
fi

log "=== Step 3 done ==="
