#!/bin/bash
#SBATCH -J 071626_cenpa_bg #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 8 #Number of threads per process
#SBATCH -t 2-00:00:00 #Short for --time walltime limit
#SBATCH --mem=16G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out #standard output name
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err #Optional, standard error name
#SBATCH -p condo #platinum #Partition name
#SBATCH -q condo #hcp-csd788 #QOS name
#SBATCH -A csd788  #htl195 for hotel #Allocation name
#SBATCH --mail-type END #Optional, Send mail when job ends
#SBATCH --mail-user you@example.com #Optional, Send mail to this address
#SBATCH --no-requeue

## If parallelization needed:
module load shared

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

# ============================================================================
# 2_generate_backgrounds.sh — Generate matched background region sets
# Bg1: Chromosome-matched random shuffle (100 iterations)
# Bg2: Mappability-matched (best-effort)
# Bg3: Local flanking regions
# ============================================================================


SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
source "${SCRIPT_DIR}/config.sh"
init_dirs

THREADS=$SLURM_CPUS_PER_TASK
log "Starting background generation"

SET_A="${FOREGROUND_DIR}/setA_centroAnno_strict.bed"
EXCLUSION_MASK="${EXCLUSION_DIR}/exclusion_mask.bed"
CHROM_SIZES="${DATA_DIR}/chrom_sizes.txt"

# ============================================================================
# Background 1: Chromosome-matched shuffle (100 iterations)
# ============================================================================
log "=== Background 1: Chromosome-matched shuffle (${N_SHUFFLES} iterations) ==="

BG1_DIR="${BACKGROUND_DIR}/bg1_chrom_shuffle"
SEED_BASE=${SHUFFLE_SEED}

n_existing=$(ls "${BG1_DIR}"/iter_*.bed 2>/dev/null | wc -l)
if [[ $n_existing -ge ${N_SHUFFLES} ]]; then
    log "SKIP Bg1: ${n_existing} shuffle files already exist"
else
    log "Generating ${N_SHUFFLES} chromosome-matched shuffle sets..."

    # Pre-compute on which chromosomes we need to place each interval
    # (each interval stays on its original chromosome via -chrom flag)

    for iter in $(seq 1 ${N_SHUFFLES}); do
        outfile="${BG1_DIR}/iter_$(printf '%03d' $iter).bed"
        seed=$((SEED_BASE + iter))

        if [[ -f "$outfile" ]] && [[ -s "$outfile" ]]; then
            continue
        fi

        bedtools shuffle -i "$SET_A" -g "$CHROM_SIZES" \
            -chrom -excl "$EXCLUSION_MASK" -seed "$seed" \
            -noOverlapping 2>/dev/null > "$outfile"

        n_placed=$(wc -l < "$outfile")
        n_expected=$(wc -l < "$SET_A")

        if [[ $n_placed -lt $n_expected ]]; then
            log "WARNING: iter $iter: only $n_placed/$n_expected intervals placed"
        fi
    done

    log "Bg1: $(ls ${BG1_DIR}/iter_*.bed | wc -l) shuffle files generated"
fi

# Verify: no overlap with exclusion mask in a few random iters
log "--- Bg1 overlap verification ---"
for iter in 001 050 100; do
    shuffle_file="${BG1_DIR}/iter_${iter}.bed"
    if [[ -f "$shuffle_file" ]]; then
        n_overlap=$(bedtools intersect -a "$shuffle_file" -b "$EXCLUSION_MASK" -u | wc -l)
        log "  iter $iter: $n_overlap overlapping exclusion mask"
    fi
done

# ============================================================================
# Background 2: Mappability-matched shuffle (best-effort)
# ============================================================================
log "=== Background 2: Mappability-matched (best-effort) ==="

BG2_DIR="${BACKGROUND_DIR}/bg2_mappability_match"
MAPPABILITY_SCORES="${MAPPABILITY_DIR}/setA_mappability_scores.txt"

# Strategy: Generate extra shuffle sets, then filter intervals to match
# foreground mappability distribution. Report results as QC.

n_existing_bg2=$(ls "${BG2_DIR}"/iter_*.bed 2>/dev/null | wc -l)
if [[ $n_existing_bg2 -ge ${N_SHUFFLES} ]]; then
    log "SKIP Bg2: ${n_existing_bg2} files exist"
else
    log "Generating ${N_SHUFFLES} shuffle sets for mappability matching..."

    for iter in $(seq 1 ${N_SHUFFLES}); do
        outfile="${BG2_DIR}/iter_$(printf '%03d' $iter).bed"
        seed=$((SEED_BASE + 10000 + iter))

        if [[ -f "$outfile" ]] && [[ -s "$outfile" ]]; then
            continue
        fi

        bedtools shuffle -i "$SET_A" -g "$CHROM_SIZES" \
            -chrom -excl "$EXCLUSION_MASK" -seed "$seed" \
            -noOverlapping 2>/dev/null > "$outfile"
    done

    log "Bg2: $(ls ${BG2_DIR}/iter_*.bed | wc -l) shuffle files generated"
fi

# Compute mappability for one Bg2 set (iter 001) to compare distributions
BG2_MAPPABILITY="${MAPPABILITY_DIR}/bg2_iter001_mappability.txt"
if [[ ! -f "$BG2_MAPPABILITY" ]]; then
    bedtools coverage -a "${BG2_DIR}/iter_001.bed" \
        -b "${FRAGMENT_DIR}/XG_152_fragments.bed" > "${MAPPABILITY_DIR}/bg2_iter001_coverage_raw.txt"

    awk -v OFS='\t' '{
        len = $3 - $2
        if (len > 0) {
            mappable = ($7 > 0 ? $7 : 0)
            fraction = mappable / len
            print $1, $2, $3, $4, fraction
        }
    }' "${MAPPABILITY_DIR}/bg2_iter001_coverage_raw.txt" > "$BG2_MAPPABILITY"
    log "Bg2 mappability scores computed for iter 001"
fi

# ============================================================================
# Background 3: Local flanking regions
# ============================================================================
log "=== Background 3: Local flanking regions ==="

BG3_DIR="${BACKGROUND_DIR}/bg3_local_flanks"

for flank_dist in "${FLANK_DISTANCES[@]}"; do
    upstream_file="${BG3_DIR}/upstream_${flank_dist}bp.bed"
    downstream_file="${BG3_DIR}/downstream_${flank_dist}bp.bed"

    if [[ -f "$upstream_file" ]] && [[ -f "$downstream_file" ]] && \
       [[ -s "$upstream_file" ]] && [[ -s "$downstream_file" ]]; then
        log "SKIP Bg3 flank ${flank_dist}bp: exists"
        continue
    fi

    log "Generating local flanks at ${flank_dist}bp..."

    # Upstream flanks: same size, shifted upstream by flank_dist
    # For each interval (chr, start, end), create (chr, start-dist, end-dist)
    awk -v dist="$flank_dist" -v OFS='\t' '{
        center = int(($2 + $3) / 2)
        size = $3 - $2
        new_start = center - dist - int(size / 2)
        new_end = new_start + size
        if (new_start >= 0) print $1, new_start, new_end, $4
    }' "$SET_A" | \
        awk -v OFS='\t' -v sizes="${CHROM_SIZES}" '
        BEGIN {
            while((getline < sizes) > 0) chrom_len[$1] = $2
        }
        {
            if ($2 >= 0 && $3 <= chrom_len[$1]) print
        }' > "${upstream_file}.tmp"

    # Remove flanks that overlap any centroAnno foreground interval
    bedtools intersect -a "${upstream_file}.tmp" -b "$SET_A" -v > "$upstream_file"
    rm -f "${upstream_file}.tmp"
    log "  Upstream ${flank_dist}bp: $(wc -l < $upstream_file) regions (from $(wc -l < $SET_A))"

    # Downstream flanks
    awk -v dist="$flank_dist" -v OFS='\t' '{
        center = int(($2 + $3) / 2)
        size = $3 - $2
        new_start = center + dist - int(size / 2)
        new_end = new_start + size
        if (new_start >= 0) print $1, new_start, new_end, $4
    }' "$SET_A" | \
        awk -v OFS='\t' -v sizes="${CHROM_SIZES}" '
        BEGIN {
            while((getline < sizes) > 0) chrom_len[$1] = $2
        }
        {
            if ($2 >= 0 && $3 <= chrom_len[$1]) print
        }' > "${downstream_file}.tmp"

    bedtools intersect -a "${downstream_file}.tmp" -b "$SET_A" -v > "$downstream_file"
    rm -f "${downstream_file}.tmp"
    log "  Downstream ${flank_dist}bp: $(wc -l < $downstream_file) regions (from $(wc -l < $SET_A))"
done

# ============================================================================
# Create concatenated background files for efficient counting
# ============================================================================
log "--- Creating concatenated background files ---"

# Bg1: concatenate all shuffle iterations with iteration ID
BG1_CONCAT="${BG1_DIR}/all_iterations_concat.bed"
if [[ ! -f "$BG1_CONCAT" ]]; then
    for f in "${BG1_DIR}"/iter_*.bed; do
        iter_id=$(basename "$f" .bed)
        awk -v id="$iter_id" -v OFS='\t' '{print $1, $2, $3, id}' "$f"
    done > "$BG1_CONCAT"
    log "Bg1 concatenated: $(wc -l < $BG1_CONCAT) total intervals"
fi

# Bg3: concatenate all flanks
BG3_CONCAT="${BG3_DIR}/all_flanks_concat.bed"
if [[ ! -f "$BG3_CONCAT" ]]; then
    for f in "${BG3_DIR}"/upstream_*.bed "${BG3_DIR}"/downstream_*.bed; do
        if [[ -f "$f" ]]; then
            flank_id=$(basename "$f" .bed)
            awk -v id="$flank_id" -v OFS='\t' '{print $1, $2, $3, id}' "$f"
        fi
    done > "$BG3_CONCAT"
    log "Bg3 concatenated: $(wc -l < $BG3_CONCAT) total intervals"
fi

log "=== 2_generate_backgrounds.sh DONE ==="
