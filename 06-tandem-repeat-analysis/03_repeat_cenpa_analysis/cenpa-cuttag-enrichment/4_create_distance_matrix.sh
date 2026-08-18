#!/bin/bash
#SBATCH -J 071826_cenpa_dist
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 1-00:00:00
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
# 4_create_distance_matrix.sh (V2) — Domain-centered metaprofile
#
# V2 changes:
#   - Center on merged domain boundaries (not individual intervals)
#   - Extended flanks: ±1 Mb
#   - Log-spaced bins relative to domain boundaries
#   - Domain interior split into margin (0-50kb from edge) and core
#   - All 4 samples including XG_153
# ============================================================================

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
source "${SCRIPT_DIR}/config.sh"
init_dirs

THREADS=$SLURM_CPUS_PER_TASK
log "Starting V2 domain-centered distance matrix generation"

MERGED_DOMAINS="${DOMAIN_DIR}/merged_domains_d${DOMAIN_MERGE_DIST}.bed"
CHROM_SIZES="${DATA_DIR}/chrom_sizes.txt"
GAP_BED="${EXCLUSION_DIR}/assembly_gaps.bed"
SET_A="${FOREGROUND_DIR}/setA_centroAnno_strict.bed"

if [[ ! -f "$MERGED_DOMAINS" ]]; then
    log "ERROR: Merged domains not found: $MERGED_DOMAINS"
    log "Run 2b_define_domains.sh first."
    exit 1
fi

mkdir -p "$DOMAIN_DISTANCE_DIR"

# ============================================================================
# Step 1: Generate boundary-relative bins for each domain
# ============================================================================
log "=== Step 1: Generating boundary-relative bins ==="

DOMAIN_BINS="${DOMAIN_DISTANCE_DIR}/all_domain_bins.bed"

if [[ -f "$DOMAIN_BINS" ]] && [[ -s "$DOMAIN_BINS" ]]; then
    log "SKIP: Domain bins already exist ($(wc -l < $DOMAIN_BINS) entries)"
else
    log "Building boundary-relative bins for $(wc -l < $MERGED_DOMAINS) domains..."

    # Read chrom sizes
    declare -A CHROM_LEN
    while IFS=$'\t' read -r chr len; do
        CHROM_LEN["$chr"]=$len
    done < "$CHROM_SIZES"

    > "$DOMAIN_BINS"

    domain_num=0
    while IFS=$'\t' read -r chr start end domain_id size; do
        domain_num=$((domain_num + 1))
        chr_max=${CHROM_LEN[$chr]:-999999999}

        # Domain midpoint and half-width
        mid=$(( (start + end) / 2 ))
        half_w=$(( (end - start) / 2 ))

        # ---- INSIDE DOMAIN (boundary-relative, negative distance) ----

        # Domain core: midpoint ± 50kb (or whole domain if < 100kb)
        core_start=$(( mid - 50000 ))
        core_end=$(( mid + 50000 ))
        [[ $core_start -lt $start ]] && core_start=$start
        [[ $core_end -gt $end ]] && core_end=$end
        if [[ $core_end -gt $core_start ]]; then
            echo -e "$chr\t$core_start\t$core_end\tdomain_${domain_num}\tdomain_core"
        fi

        # Left margin: 0-50kb from left boundary (inside)
        left_margin_end=$(( start + 50000 ))
        [[ $left_margin_end -gt $core_start ]] && left_margin_end=$core_start
        if [[ $left_margin_end -gt $start ]]; then
            echo -e "$chr\t$start\t$left_margin_end\tdomain_${domain_num}\tleft_margin_0_50kb"
        fi

        # Right margin: 0-50kb from right boundary (inside)
        right_margin_start=$(( end - 50000 ))
        [[ $right_margin_start -lt $core_end ]] && right_margin_start=$core_end
        if [[ $end -gt $right_margin_start ]]; then
            echo -e "$chr\t$right_margin_start\t$end\tdomain_${domain_num}\tright_margin_0_50kb"
        fi

        # ---- OUTSIDE DOMAIN: LEFT FLANKS (upstream) ----
        # left_0_1kb:   [start-1000, start]
        # left_1_10kb:  [start-10000, start-1000]
        # left_10_100kb:[start-100000, start-10000]
        # left_100_500kb:[start-500000, start-100000]
        # left_500kb_1Mb:[start-1000000, start-500000]

        for bin_def in "0_1000:1000" "1000_10000:10000" "10000_100000:100000" "100000_500000:500000" "500000_1000000:1000000"; do
            bin_label="left_${bin_def%:*}"
            bin_label="${bin_label//_/-}"  # replace first _
            bin_dist=${bin_def#*:}

            inner_dist=${bin_def%:*}
            outer_dist=${inner_dist#*_}
            inner_dist=${inner_dist%_*}

            b_start=$(( start - outer_dist ))
            b_end=$(( start - inner_dist ))
            [[ $b_start -lt 0 ]] && b_start=0
            if [[ $b_end -gt $b_start ]]; then
                echo -e "$chr\t$b_start\t$b_end\tdomain_${domain_num}\tleft_${inner_dist}_${outer_dist}"
            fi
        done

        # ---- OUTSIDE DOMAIN: RIGHT FLANKS (downstream) ----
        for bin_def in "0_1000:1000" "1000_10000:10000" "10000_100000:100000" "100000_500000:500000" "500000_1000000:1000000"; do
            bin_label="right_${bin_def%:*}"
            bin_dist=${bin_def#*:}

            inner_dist=${bin_def%:*}
            outer_dist=${inner_dist#*_}
            inner_dist=${inner_dist%_*}

            b_start=$(( end + inner_dist ))
            b_end=$(( end + outer_dist ))
            [[ $b_end -gt $chr_max ]] && b_end=$chr_max
            if [[ $b_end -gt $b_start ]]; then
                echo -e "$chr\t$b_start\t$b_end\tdomain_${domain_num}\tright_${inner_dist}_${outer_dist}"
            fi
        done

    done < "$MERGED_DOMAINS" >> "$DOMAIN_BINS"

    n_bins=$(wc -l < "$DOMAIN_BINS")
    log "Generated $n_bins boundary-relative bins for $domain_num domains"

    # Per-bin-type counts
    log "Bin type counts:"
    awk '{print $5}' "$DOMAIN_BINS" | sort | uniq -c | sort -rn | while read count label; do
        log "  $label: $count"
    done
fi

# ============================================================================
# Step 2: Count fragments in each bin
# ============================================================================
log "=== Step 2: Counting fragments in domain bins ==="

# Get chromosomes with domains
DOMAIN_CHROMS=$(cut -f1 "$MERGED_DOMAINS" | sort -u | tr '\n' ' ')
log "Chromosomes with domains: $DOMAIN_CHROMS"

for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    count_file="${DOMAIN_DISTANCE_DIR}/${sample}_domain_bin_counts.txt"

    if [[ -f "$count_file" ]] && [[ -s "$count_file" ]]; then
        log "SKIP: $sample domain bin counts exist ($(wc -l < $count_file) lines)"
        continue
    fi

    # Subset fragment BED to domain chromosomes
    frag_bed_subset="${DOMAIN_DISTANCE_DIR}/${sample}_fragments_domainChroms.bed"
    if [[ ! -f "$frag_bed_subset" ]] || [[ ! -s "$frag_bed_subset" ]]; then
        log "Subsetting $sample fragment BED to domain chromosomes..."
        awk -v chrs="$DOMAIN_CHROMS" \
            'BEGIN{split(chrs, a, " "); for(i in a) keep[a[i]]=1}
             $1 in keep' \
            "$frag_bed" > "$frag_bed_subset" || true
        log "  Reduced: $(wc -l < $frag_bed_subset) / $(wc -l < $frag_bed) lines"
    fi

    log "Counting $sample fragments in domain bins..."
    bedtools coverage -a "$DOMAIN_BINS" -b "$frag_bed_subset" -counts > "$count_file"
    log "  Done: $(wc -l < $count_file) lines"
done

# ============================================================================
# Step 3: Also generate V1-style interval bins (for backward compatibility)
# ============================================================================
log "=== Step 3: V1-style interval bins (for comparison) ==="

DISTANCE_DIR="${DATA_DIR}/distance_profiles"
mkdir -p "$DISTANCE_DIR"

FLANK_BINS_ALL="${DISTANCE_DIR}/all_flank_bins.bed"

if [[ -f "$FLANK_BINS_ALL" ]] && [[ -s "$FLANK_BINS_ALL" ]]; then
    log "SKIP: V1 flank bins already exist ($(wc -l < $FLANK_BINS_ALL) entries)"
else
    log "Building V1 flank bins (interval-centered, XG_153 now included)..."

    declare -A CHROM_LEN_V1
    while IFS=$'\t' read -r chr len; do
        CHROM_LEN_V1["$chr"]=$len
    done < "$CHROM_SIZES"

    BINS=("-1000000:-500000" "-500000:-100000" "-100000:-10000" "-10000:-1000" "-1000:0" "body" "0:1000" "1000:10000" "10000:100000" "100000:500000" "500000:1000000")
    BIN_LABELS=("up_1M_500k" "up_500k_100k" "up_100k_10k" "up_10k_1k" "up_1k_0" "body" "down_0_1k" "down_1k_10k" "down_10k_100k" "down_100k_500k" "down_500k_1M")

    > "$FLANK_BINS_ALL"
    interval_id=0

    while IFS=$'\t' read -r chr start end orig_id; do
        interval_id=$((interval_id + 1))
        for i in "${!BINS[@]}"; do
            bin="${BINS[$i]}"
            label="${BIN_LABELS[$i]}"

            if [[ "$bin" == "body" ]]; then
                echo -e "$chr\t$start\t$end\tinterval_${interval_id}\t${label}"
                continue
            fi

            bin_start=${bin%:*}
            bin_end=${bin#*:}

            if [[ $bin_start -lt 0 ]]; then
                b_start=$((start + bin_start))
                b_end=$((start + bin_end))
                [[ $b_start -lt 0 ]] && b_start=0
            elif [[ $bin_end -gt 0 ]]; then
                b_start=$((end + bin_start))
                b_end=$((end + bin_end))
                [[ $b_end -gt ${CHROM_LEN_V1[$chr]:-999999999} ]] && b_end=${CHROM_LEN_V1[$chr]}
            fi

            if [[ ${b_end:-0} -gt ${b_start:-0} ]]; then
                echo -e "$chr\t$b_start\t$b_end\tinterval_${interval_id}\t${label}"
            fi
        done
    done < "$SET_A" >> "$FLANK_BINS_ALL"

    log "Generated $(wc -l < $FLANK_BINS_ALL) V1 flank bins"
fi

# Count fragments in V1 bins (for all samples including XG_153)
CENTRO_CHROMS=$(cut -f1 "$SET_A" | sort -u | tr '\n' ' ')
for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    count_file="${DISTANCE_DIR}/${sample}_distance_bins_counts.txt"

    if [[ -f "$count_file" ]] && [[ -s "$count_file" ]]; then
        log "SKIP: $sample V1 distance counts exist ($(wc -l < $count_file) lines)"
        continue
    fi

    frag_bed_subset="${DISTANCE_DIR}/${sample}_fragments_centroChroms.bed"
    if [[ ! -f "$frag_bed_subset" ]] || [[ ! -s "$frag_bed_subset" ]]; then
        awk -v chrs="$CENTRO_CHROMS" \
            'BEGIN{split(chrs, a, " "); for(i in a) keep[a[i]]=1}
             $1 in keep' \
            "$frag_bed" > "$frag_bed_subset" || true
    fi

    log "Counting $sample in V1 distance bins..."
    bedtools coverage -a "$FLANK_BINS_ALL" -b "$frag_bed_subset" -counts > "$count_file"
    log "  Done: $(wc -l < $count_file) lines"
done

log "=== 4_create_distance_matrix.sh DONE ==="
