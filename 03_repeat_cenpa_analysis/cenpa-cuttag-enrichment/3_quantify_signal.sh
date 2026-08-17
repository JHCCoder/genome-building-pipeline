#!/bin/bash
#SBATCH -J 071826_cenpa_quant
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 2-00:00:00
#SBATCH --mem=64G
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user you@example.com
#SBATCH --no-requeue

## If parallelization needed:
module load shared

source /tscc/nfs/home/jhc103/.bashrc
conda activate bulk-HiC-processing

# ============================================================================
# 3_quantify_signal.sh (V2) — bedtools coverage on all region sets × all samples
#
# V2 additions:
#   - Merged domain foreground counting
#   - Domain local flank background (±500 kb, ±1 Mb)
#   - XG_153 integrated throughout
#   - Per-base coverage tracks for domain-centered metaprofiles
# ============================================================================

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
source "${SCRIPT_DIR}/config.sh"
init_dirs

THREADS=$SLURM_CPUS_PER_TASK
log "Starting V2 quantification (threads=$THREADS)"

SET_A="${FOREGROUND_DIR}/setA_centroAnno_strict.bed"
CHROM_SIZES="${DATA_DIR}/chrom_sizes.txt"
GAP_BED="${EXCLUSION_DIR}/assembly_gaps.bed"
MERGED_DOMAINS="${DOMAIN_DIR}/merged_domains_d${DOMAIN_MERGE_DIST}.bed"

# Check that merged domains exist
if [[ ! -f "$MERGED_DOMAINS" ]]; then
    log "ERROR: Merged domains not found: $MERGED_DOMAINS"
    log "Run 2b_define_domains.sh first."
    exit 1
fi
log "Merged domains: $(wc -l < $MERGED_DOMAINS) domains"

# ============================================================================
# Define all region sets to quantify
# ============================================================================
declare -A REGION_SETS
REGION_SETS["setA_strict"]="$SET_A"
REGION_SETS["setB_d50k_p10k"]="${FOREGROUND_DIR}/setB_merged_d${MERGE_DIST_1}_pad${FLANK_PAD_1}.bed"
REGION_SETS["setB_d25k_p5k"]="${FOREGROUND_DIR}/setB_merged_d${MERGE_DIST_2}_pad${FLANK_PAD_2}.bed"
REGION_SETS["setB_d100k_p20k"]="${FOREGROUND_DIR}/setB_merged_d${MERGE_DIST_3}_pad${FLANK_PAD_3}.bed"
REGION_SETS["setC_hicat_full"]="${FOREGROUND_DIR}/setC_hicat_chr4_full.bed"
REGION_SETS["setC_hicat_occupied"]="${FOREGROUND_DIR}/setC_hicat_chr4_occupied.bed"
REGION_SETS["setA_chr4"]="${FOREGROUND_DIR}/setA_chr4_only.bed"

# Background concatenated files
BG1_CONCAT="${BACKGROUND_DIR}/bg1_chrom_shuffle/all_iterations_concat.bed"
BG3_CONCAT="${BACKGROUND_DIR}/bg3_local_flanks/all_flanks_concat.bed"

# ============================================================================
# SECTION 1: Domain local flank generation
# ============================================================================
log "=== SECTION 1: Domain local flank generation ==="

for flank_dist in 500000 1000000; do
    flank_label=$(echo $flank_dist | awk '{printf "%.0fkb", $1/1000}')
    upstream_file="${DOMAIN_FLANKS_DIR}/upstream_${flank_label}.bed"
    downstream_file="${DOMAIN_FLANKS_DIR}/downstream_${flank_label}.bed"

    if [[ -f "$upstream_file" ]] && [[ -f "$downstream_file" ]] && \
       [[ -s "$upstream_file" ]] && [[ -s "$downstream_file" ]]; then
        log "SKIP domain flanks ${flank_label}: exist"
        continue
    fi

    log "Generating domain flanks at ${flank_label}..."

    # Upstream: shift each domain upstream by flank_dist, same size
    awk -v dist="$flank_dist" -v OFS='\t' '{
        center = int(($2 + $3) / 2)
        size = $3 - $2
        new_start = center - dist - int(size / 2)
        new_end = new_start + size
        if (new_start >= 0) print $1, new_start, new_end, $4
    }' "$MERGED_DOMAINS" | \
        awk -v OFS='\t' -v sizes="${CHROM_SIZES}" '
        BEGIN {
            while((getline < sizes) > 0) chrom_len[$1] = $2
        }
        {
            if ($2 >= 0 && $3 <= chrom_len[$1]) print
        }' > "${upstream_file}.tmp"

    # Exclude overlaps with centroAnno intervals
    bedtools intersect -a "${upstream_file}.tmp" -b "$SET_A" -v > "${upstream_file}.tmp2"
    # Exclude gaps
    bedtools intersect -a "${upstream_file}.tmp2" -b "$GAP_BED" -v > "$upstream_file"
    rm -f "${upstream_file}.tmp" "${upstream_file}.tmp2"
    log "  Upstream ${flank_label}: $(wc -l < $upstream_file) regions"

    # Downstream: shift downstream by flank_dist
    awk -v dist="$flank_dist" -v OFS='\t' '{
        center = int(($2 + $3) / 2)
        size = $3 - $2
        new_start = center + dist - int(size / 2)
        new_end = new_start + size
        if (new_start >= 0) print $1, new_start, new_end, $4
    }' "$MERGED_DOMAINS" | \
        awk -v OFS='\t' -v sizes="${CHROM_SIZES}" '
        BEGIN {
            while((getline < sizes) > 0) chrom_len[$1] = $2
        }
        {
            if ($2 >= 0 && $3 <= chrom_len[$1]) print
        }' > "${downstream_file}.tmp"

    bedtools intersect -a "${downstream_file}.tmp" -b "$SET_A" -v > "${downstream_file}.tmp2"
    bedtools intersect -a "${downstream_file}.tmp2" -b "$GAP_BED" -v > "$downstream_file"
    rm -f "${downstream_file}.tmp" "${downstream_file}.tmp2"
    log "  Downstream ${flank_label}: $(wc -l < $downstream_file) regions"
done

# Concatenate all domain flanks for efficient counting
DOMAIN_FLANKS_CONCAT="${DOMAIN_FLANKS_DIR}/all_domain_flanks.bed"
if [[ ! -f "$DOMAIN_FLANKS_CONCAT" ]]; then
    for f in "${DOMAIN_FLANKS_DIR}"/upstream_*.bed "${DOMAIN_FLANKS_DIR}"/downstream_*.bed; do
        if [[ -f "$f" ]]; then
            flank_id=$(basename "$f" .bed)
            awk -v id="$flank_id" -v OFS='\t' '{print $1, $2, $3, id}' "$f"
        fi
    done > "$DOMAIN_FLANKS_CONCAT"
    log "Domain flanks concatenated: $(wc -l < $DOMAIN_FLANKS_CONCAT) total regions"
fi

# ============================================================================
# SECTION 2: Fragment counting — existing foreground sets (V1 + XG_153)
# ============================================================================
log "=== SECTION 2: Foreground fragment counting (Set A/B/C) ==="

for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    check_file "$frag_bed" "Fragment BED $sample"

    for rs_name in "${!REGION_SETS[@]}"; do
        rs_file="${REGION_SETS[$rs_name]}"
        count_file="${COUNTS_DIR}/${sample}_${rs_name}.txt"

        if [[ -f "$count_file" ]] && [[ -s "$count_file" ]]; then
            log "SKIP: ${sample}_${rs_name} ($(wc -l < $count_file) lines)"
            continue
        fi

        if [[ ! -f "$rs_file" ]] || [[ ! -s "$rs_file" ]]; then
            log "SKIP: ${sample}_${rs_name} (region file empty/missing)"
            continue
        fi

        log "Counting: $sample x $rs_name ($(wc -l < $rs_file) regions)"

        bedtools coverage -a "$rs_file" -b "$frag_bed" -counts > "$count_file"
        total_frags=$(awk '{sum+=$4}END{print sum}' "$count_file")
        log "  Total fragments in regions: $total_frags"
    done
done

# ============================================================================
# SECTION 3: Fragment counting — Merged domains (V2 foreground)
# ============================================================================
log "=== SECTION 3: Domain foreground counting ==="

for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    count_file="${DOMAIN_COUNTS_DIR}/${sample}_domains.txt"

    if [[ -f "$count_file" ]] && [[ -s "$count_file" ]]; then
        log "SKIP: ${sample}_domains ($(wc -l < $count_file) lines)"
        continue
    fi

    log "Counting: $sample x merged domains ($(wc -l < $MERGED_DOMAINS) domains)"
    bedtools coverage -a "$MERGED_DOMAINS" -b "$frag_bed" -counts > "$count_file"
    log "  Done: $(wc -l < $count_file) domains"
done

# ============================================================================
# SECTION 4: Fragment counting — Domain local flanks (V2 local background)
# ============================================================================
log "=== SECTION 4: Domain flank counting (local background) ==="

for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    count_file="${DOMAIN_COUNTS_DIR}/${sample}_domain_flanks.txt"

    if [[ -f "$count_file" ]] && [[ -s "$count_file" ]]; then
        log "SKIP: ${sample}_domain_flanks ($(wc -l < $count_file) lines)"
        continue
    fi

    if [[ ! -f "$DOMAIN_FLANKS_CONCAT" ]]; then
        log "WARNING: Domain flanks concatenated file not found, counting individually..."
        for flank_file in "${DOMAIN_FLANKS_DIR}"/upstream_*.bed \
                          "${DOMAIN_FLANKS_DIR}"/downstream_*.bed; do
            if [[ -f "$flank_file" ]]; then
                flank_id=$(basename "$flank_file" .bed)
                flank_count="${DOMAIN_COUNTS_DIR}/${sample}_domain_flank_${flank_id}.txt"
                if [[ ! -f "$flank_count" ]]; then
                    bedtools coverage -a "$flank_file" -b "$frag_bed" -counts > "$flank_count"
                fi
            fi
        done
        continue
    fi

    log "Counting: $sample x domain flanks ($(wc -l < $DOMAIN_FLANKS_CONCAT) total)"
    bedtools coverage -a "$DOMAIN_FLANKS_CONCAT" -b "$frag_bed" -counts > "$count_file"
    log "  Done: $(wc -l < $count_file) lines"
done

# ============================================================================
# SECTION 5: Fragment counting — Bg1 (chromosome shuffle, all iterations)
# ============================================================================
log "=== SECTION 5: Background 1 counting (shuffled) ==="

for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    count_file="${COUNTS_DIR}/${sample}_bg1_chrom_shuffle.txt"

    if [[ -f "$count_file" ]] && [[ -s "$count_file" ]]; then
        log "SKIP: ${sample}_bg1 ($(wc -l < $count_file) lines)"
        continue
    fi

    if [[ ! -f "$BG1_CONCAT" ]]; then
        log "WARNING: Bg1 concatenated file not found, counting individually..."
        for iter_file in "${BACKGROUND_DIR}/bg1_chrom_shuffle"/iter_*.bed; do
            iter_id=$(basename "$iter_file" .bed)
            iter_count="${COUNTS_DIR}/${sample}_bg1_${iter_id}.txt"
            if [[ ! -f "$iter_count" ]]; then
                bedtools coverage -a "$iter_file" -b "$frag_bed" -counts > "$iter_count"
            fi
        done
        continue
    fi

    log "Counting: $sample x Bg1 ($(wc -l < $BG1_CONCAT) total intervals)"
    bedtools coverage -a "$BG1_CONCAT" -b "$frag_bed" -counts > "$count_file"
    log "  Done: $(wc -l < $count_file) lines"
done

# ============================================================================
# SECTION 6: Fragment counting — Bg3 (local flanks, V1)
# ============================================================================
log "=== SECTION 6: Background 3 counting (local flanks) ==="

for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    count_file="${COUNTS_DIR}/${sample}_bg3_local_flanks.txt"

    if [[ -f "$count_file" ]] && [[ -s "$count_file" ]]; then
        log "SKIP: ${sample}_bg3 ($(wc -l < $count_file) lines)"
        continue
    fi

    if [[ -f "$BG3_CONCAT" ]]; then
        log "Counting: $sample x Bg3 ($(wc -l < $BG3_CONCAT) total intervals)"
        bedtools coverage -a "$BG3_CONCAT" -b "$frag_bed" -counts > "$count_file"
        log "  Done: $(wc -l < $count_file) lines"
    fi
done

# ============================================================================
# SECTION 7: Per-base coverage tracks (for domain metaprofiles)
# ============================================================================
log "=== SECTION 7: Per-base coverage for domain metaprofiles ==="

for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    cov_file="${COVERAGE_DIR}/${sample}_perbase.txt.gz"

    if [[ -f "$cov_file" ]] && [[ -s "$cov_file" ]]; then
        log "SKIP per-base coverage $sample: exists"
        continue
    fi

    log "Generating per-base coverage for $sample..."
    bedtools genomecov -i "$frag_bed" -g "$CHROM_SIZES" -d | gzip > "$cov_file"
    log "  Done: $(zcat $cov_file | wc -l) positions"
done

# ============================================================================
# Summary
# ============================================================================
log "=== Quantification summary (V2) ==="
log "Foreground count files:"
ls -la "${COUNTS_DIR}"/*.txt 2>/dev/null | awk '{print "  " $NF " (" $5 " bytes)"}' | tail -20
log ""
log "Domain count files:"
ls -la "${DOMAIN_COUNTS_DIR}"/*.txt 2>/dev/null | awk '{print "  " $NF " (" $5 " bytes)"}' | tail -20
log ""
log "Coverage files:"
ls -la "${COVERAGE_DIR}"/*.txt.gz 2>/dev/null | awk '{print "  " $NF " (" $5 " bytes)"}'
log ""
log "Domain flanks:"
ls -la "${DOMAIN_FLANKS_DIR}"/*.bed 2>/dev/null | awk '{print "  " $NF " (" $5 " bytes)"}' | tail -10

log "=== 3_quantify_signal.sh DONE ==="
