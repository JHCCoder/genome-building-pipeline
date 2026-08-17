#!/bin/bash
#SBATCH -J 071626_cenpa_prep #Optional, short for --job-name
#SBATCH -N 1 #Number of nodes
#SBATCH -n 1 #Total number of tasks increase this number to increase parallelization
#SBATCH -c 16 #Number of threads per process
#SBATCH -t 1-00:00:00 #Short for --time walltime limit
#SBATCH --mem=32G
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

SCRIPT_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
source "${SCRIPT_DIR}/config.sh"
init_dirs

THREADS=$SLURM_CPUS_PER_TASK
log "Starting validation and preparation (threads=$THREADS)"

# ============================================================================
# SECTION 1: Input validation
# ============================================================================
log "=== SECTION 1: Input validation ==="

# Check BAMs and indices
for sample in "${BAM_LIST[@]}"; do
    bam="${BAM_PATHS[$sample]}"
    check_file "$bam" "BAM $sample"
    check_file "${bam}.bai" "BAM index $sample"
done

# Check centroAnno files
check_file "$REPEAT_REGIONS" "centroAnno repeat_regions.bed"
check_file "$HORS_BED" "centroAnno HORs.bed"
check_file "$TOP10_REPEATS" "centroAnno top10_repeats.bed"

# Check assembly
check_file "$ASSEMBLY_FASTA" "Assembly FASTA"
check_file "$ASSEMBLY_FAI" "Assembly FAI"

# Check HiCAT chr4
check_file "$HICAT_DECOMP" "HiCAT final_decomposition.bed"

log "All input files present"

# ============================================================================
# SECTION 1b: Read counts and flag validation
# ============================================================================
log "--- Read counts per sample per chromosome ---"
for sample in "${BAM_LIST[@]}"; do
    bam="${BAM_PATHS[$sample]}"
    log "Sample $sample:"
    samtools idxstats "$bam" 2>/dev/null | head -35 > "${QC_DIR}/${sample}_idxstats.txt"
    total_reads=$(awk '{sum+=$3} END{print sum}' "${QC_DIR}/${sample}_idxstats.txt")
    log "  Total mapped reads: $total_reads"
    echo "$sample $total_reads" >> "${QC_DIR}/library_sizes_reads.txt"
done

log "--- Flag distribution checks ---"
for sample in "${BAM_LIST[@]}"; do
    bam="${BAM_PATHS[$sample]}"
    flags=$(samtools view "$bam" 2>/dev/null | head -100000 | awk '{print $2}' | sort | uniq -c | sort -rn)
    echo "=== $sample flags ===" > "${QC_DIR}/${sample}_flags.txt"
    echo "$flags" >> "${QC_DIR}/${sample}_flags.txt"
    log "$sample flags: $(echo "$flags" | head -4 | tr '\n' ' ')"

    # Check for unexpected flags
    unexpected=$(echo "$flags" | awk '$2!=99 && $2!=147 && $2!=163 && $2!=83 {print $2}')
    if [[ -n "$unexpected" ]]; then
        log "WARNING: Unexpected flags in $sample: $unexpected"
    else
        log "OK: $sample only has expected flags (99/147/163/83)"
    fi

    # Check for secondary (256) or supplementary (2048)
    secondary_count=$(samtools view -c -f 0x100 "$bam" 2>/dev/null)
    supplementary_count=$(samtools view -c -f 0x800 "$bam" 2>/dev/null)
    log "  Secondary alignments: $secondary_count, Supplementary: $supplementary_count"
done

# ============================================================================
# SECTION 1c: Chromosome naming consistency check
# ============================================================================
log "--- Chromosome naming consistency ---"

# Extract chromosome names from BAM
samtools idxstats "${BAM_PATHS["XG_150"]}" 2>/dev/null | awk '$1 ~ /^chr/{print $1}' | sort > "${QC_DIR}/bam_chroms.txt"

# Extract chromosome names from centroAnno BED
awk '{print $1}' "$REPEAT_REGIONS" | sort -u > "${QC_DIR}/centroAnno_chroms.txt"

# Extract chromosome names from FAI
awk '{print $1}' "$ASSEMBLY_FAI" | grep '^chr' | sort > "${QC_DIR}/fasta_chroms.txt"

# Check for mismatches
comm -23 "${QC_DIR}/centroAnno_chroms.txt" "${QC_DIR}/bam_chroms.txt" > "${QC_DIR}/chrom_mismatch_centroAnno_vs_bam.txt"
comm -23 "${QC_DIR}/bam_chroms.txt" "${QC_DIR}/centroAnno_chroms.txt" > "${QC_DIR}/chrom_mismatch_bam_vs_centroAnno.txt"

if [[ -s "${QC_DIR}/chrom_mismatch_centroAnno_vs_bam.txt" ]]; then
    log "WARNING: centroAnno chroms not in BAM:"
    cat "${QC_DIR}/chrom_mismatch_centroAnno_vs_bam.txt"
fi
if [[ -s "${QC_DIR}/chrom_mismatch_bam_vs_centroAnno.txt" ]]; then
    log "INFO: BAM chroms not in centroAnno:"
    cat "${QC_DIR}/chrom_mismatch_bam_vs_centroAnno.txt"
fi

# Verify all expected chromosomes present
for chr in $CHROMOSOMES; do
    if ! grep -qx "$chr" "${QC_DIR}/bam_chroms.txt"; then
        log "WARNING: $chr not found in BAM"
    fi
    if ! grep -qx "$chr" "${QC_DIR}/centroAnno_chroms.txt"; then
        log "WARNING: $chr not found in centroAnno BED"
    fi
done
log "OK: Chromosome naming consistency checked"

# ============================================================================
# SECTION 2: Fragment BED generation
# ============================================================================
log "=== SECTION 2: Fragment BED generation ==="

for sample in "${BAM_LIST[@]}"; do
    bam="${BAM_PATHS[$sample]}"
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"

    if [[ -f "$frag_bed" ]] && [[ -s "$frag_bed" ]]; then
        n_frags=$(wc -l < "$frag_bed")
        log "SKIP $sample: fragment BED exists ($n_frags fragments)"
        echo "$sample $n_frags" >> "${QC_DIR}/library_sizes_fragments.txt"
        continue
    fi

    log "Generating fragment BED for $sample..."

    # Filter to proper pairs, exclude secondary/supplementary, MAPQ threshold
    # Name-sort, convert to BEDPE, extract fragment coordinates
    # Output: one row per fragment (leftmost to rightmost coordinate)
    samtools view -b ${SAMTOOLS_FILTER} -q ${MAPQ_THRESHOLD} "$bam" | \
        samtools sort -n -@ $THREADS - | \
        bedtools bamtobed -bedpe -i - | \
        awk 'BEGIN{OFS="\t"}{
            if($1==$4 && $6>$2) print $1, $2, $6
        }' | \
        sort -k1,1 -k2,2n --parallel=$THREADS -S 4G > "$frag_bed"

    n_frags=$(wc -l < "$frag_bed")
    log "  Generated $n_frags fragments for $sample"
    echo "$sample $n_frags" >> "${QC_DIR}/library_sizes_fragments.txt"

    # Verify: fragment count should be ~half of read count
    n_reads=$(grep "^$sample " "${QC_DIR}/library_sizes_reads.txt" | awk '{print $2}')
    ratio=$(echo "scale=2; $n_frags / $n_reads" | bc)
    log "  Fragment/read ratio: $ratio (expected ~0.5 for paired-end)"
done

log "Fragment BEDs complete"

# ============================================================================
# SECTION 2b: Fragment size distributions (QC)
# ============================================================================
log "--- Fragment size distributions ---"
for sample in "${BAM_LIST[@]}"; do
    frag_bed="${FRAGMENT_DIR}/${sample}_fragments.bed"
    awk '{len=$3-$2; if(len>0 && len<2000) print len}' "$frag_bed" | \
        sort -n | uniq -c > "${QC_DIR}/${sample}_fragment_sizes.txt"
    median_size=$(awk '{l=$3-$2; sum+=l; count++; arr[count]=l} END{
        asort(arr); print arr[int(count/2)]
    }' "$frag_bed")
    log "  $sample median fragment size: ${median_size}bp"
done

# ============================================================================
# SECTION 3: Chromosome sizes file
# ============================================================================
log "=== SECTION 3: Chromosome sizes ==="
grep -E '^chr([0-9]+|[XY])\s' "$ASSEMBLY_FAI" | \
    awk '{print $1"\t"$2}' | sort -k1,1V > "${DATA_DIR}/chrom_sizes.txt"

# Also create chr-only sizes (for bedtools shuffle -chrom)
awk '{print $1"\t"$2}' "${DATA_DIR}/chrom_sizes.txt" > "${DATA_DIR}/chrom_sizes_chrOnly.txt"

n_chroms=$(wc -l < "${DATA_DIR}/chrom_sizes.txt")
log "Created chrom_sizes.txt with $n_chroms chromosomes"
log "Total assembly size (chr1-28,X,Y): $(awk '{sum+=$2}END{print sum}' ${DATA_DIR}/chrom_sizes.txt) bp"

# ============================================================================
# SECTION 4: Foreground Set A — centroAnno strict
# ============================================================================
log "=== SECTION 4: Foreground Set A (strict centroAnno) ==="

# Filter to chr1-28,X,Y only, sort by coordinate
awk -v OFS='\t' '{
    if($1 ~ /^chr([0-9]+|[XY])$/) {
        print $1, $2, $3, NR
    }
}' "$REPEAT_REGIONS" | sort -k1,1V -k2,2n > "${FOREGROUND_DIR}/setA_centroAnno_strict.bed"

n_setA=$(wc -l < "${FOREGROUND_DIR}/setA_centroAnno_strict.bed")
log "Set A (strict): $n_setA intervals"

# Report per-chromosome counts
awk '{print $1}' "${FOREGROUND_DIR}/setA_centroAnno_strict.bed" | sort | uniq -c | \
    sort -k1 -n > "${QC_DIR}/setA_per_chromosome_counts.txt"
log "Per-chromosome Set A counts:"
cat "${QC_DIR}/setA_per_chromosome_counts.txt"

# Check for overlapping intervals on same chromosome
log "--- Overlapping interval check ---"
n_overlaps=$(bedtools intersect -a "${FOREGROUND_DIR}/setA_centroAnno_strict.bed" \
    -b "${FOREGROUND_DIR}/setA_centroAnno_strict.bed" -wao | \
    awk '$4!=$8 && $8!="."' | wc -l)
log "  Number of overlapping interval pairs: $n_overlaps"

# ============================================================================
# SECTION 5: Foreground Set B — centroAnno repeat neighborhoods
# ============================================================================
log "=== SECTION 5: Foreground Set B (repeat neighborhoods) ==="

for variant in "1_${MERGE_DIST_1}_${FLANK_PAD_1}" "2_${MERGE_DIST_2}_${FLANK_PAD_2}" "3_${MERGE_DIST_3}_${FLANK_PAD_3}"; do
    merge_d=${variant#*_}; merge_d=${merge_d%_*}
    flank_d=${variant##*_}

    outfile="${FOREGROUND_DIR}/setB_merged_d${merge_d}_pad${flank_d}.bed"

    if [[ -f "$outfile" ]] && [[ -s "$outfile" ]]; then
        log "SKIP Set B variant d=${merge_d} pad=${flank_d}: exists ($(wc -l < $outfile) regions)"
        continue
    fi

    log "Set B variant: merge -d ${merge_d} | slop -b ${flank_d}"

    bedtools merge -d "$merge_d" -i "${FOREGROUND_DIR}/setA_centroAnno_strict.bed" | \
        bedtools slop -b "$flank_d" -g "${DATA_DIR}/chrom_sizes.txt" | \
        awk -v OFS='\t' '{if($2>=0) print $1, $2, $3, NR}' > "$outfile"

    n_setB=$(wc -l < "$outfile")
    log "  Regions after merge+pad: $n_setB"
done

# ============================================================================
# SECTION 6: Foreground Set C — chr4 HiCAT domains
# ============================================================================
log "=== SECTION 6: Foreground Set C (chr4 HiCAT) ==="

# HiCAT spans chr4:124,000,000-134,000,000 (10Mb)
# Extract this as a single HOR domain; also create sub-domains for within-HiCAT analysis

# Full HiCAT domain
echo -e "chr4\t124000000\t134000000\tHiCAT_full_domain" > "${FOREGROUND_DIR}/setC_hicat_chr4_full.bed"

# Also merge HiCAT monomers to get exact occupied regions (vs the 10Mb bounding box)
awk '$1=="chr4"{print $1"\t"$2"\t"$3}' "$HICAT_DECOMP" | \
    sort -k1,1 -k2,2n | \
    bedtools merge -d 1000 -i - > "${FOREGROUND_DIR}/setC_hicat_chr4_occupied.bed"

n_hicat_regions=$(wc -l < "${FOREGROUND_DIR}/setC_hicat_chr4_occupied.bed")
log "Set C (HiCAT chr4): $n_hicat_regions occupied regions within 124-134Mb domain"

# Also extract centroAnno chr4 seeds for direct comparison
awk '$1=="chr4"' "${FOREGROUND_DIR}/setA_centroAnno_strict.bed" > "${FOREGROUND_DIR}/setA_chr4_only.bed"
log "Set A chr4 seeds: $(wc -l < ${FOREGROUND_DIR}/setA_chr4_only.bed) intervals"

# ============================================================================
# SECTION 7: Gap BED generation
# ============================================================================
log "=== SECTION 7: Gap BED generation ==="

GAP_BED="${EXCLUSION_DIR}/assembly_gaps.bed"

if [[ -f "$GAP_BED" ]] && [[ -s "$GAP_BED" ]]; then
    log "SKIP gap BED: exists ($(wc -l < $GAP_BED) gap regions)"
else
    log "Extracting n-stretches from assembly FASTA..."
    python3 -c "
import sys

with open('${ASSEMBLY_FASTA}') as f:
    chrom = None
    pos = 0
    in_gap = False
    gap_start = 0
    for line in f:
        if line.startswith('>'):
            if in_gap:
                print(f'{chrom}\t{gap_start}\t{pos}')
                in_gap = False
            chrom = line[1:].split()[0]
            pos = 0
        else:
            for base in line.strip():
                if base.lower() == 'n':
                    if not in_gap:
                        gap_start = pos
                        in_gap = True
                else:
                    if in_gap:
                        print(f'{chrom}\t{gap_start}\t{pos}')
                        in_gap = False
                pos += 1
    if in_gap:
        print(f'{chrom}\t{gap_start}\t{pos}')
" > "$GAP_BED"

    n_gaps=$(wc -l < "$GAP_BED")
    log "Found $n_gaps gap regions"
fi

# Count gaps per chromosome
awk '{print $1}' "$GAP_BED" | sort | uniq -c | sort -k2,2V > "${QC_DIR}/gaps_per_chromosome.txt"
log "Gaps per chromosome (first 10):"
head -10 "${QC_DIR}/gaps_per_chromosome.txt"

# ============================================================================
# SECTION 8: Exclusion mask
# ============================================================================
log "=== SECTION 8: Exclusion mask ==="

EXCLUSION_MASK="${EXCLUSION_DIR}/exclusion_mask.bed"

# Combine: all foreground intervals + gaps + unplaced scaffolds
cat "${FOREGROUND_DIR}/setA_centroAnno_strict.bed" \
    "${FOREGROUND_DIR}/setB_merged_d${MERGE_DIST_1}_pad${FLANK_PAD_1}.bed" \
    "$GAP_BED" | \
    awk -v OFS='\t' '{print $1, $2, $3}' | \
    sort -k1,1V -k2,2n | \
    bedtools merge -i - > "${EXCLUSION_DIR}/foreground_and_gaps.bed"

# Add unplaced scaffolds
awk '$1 ~ /^seq/{print $1"\t0\t"$2}' "$ASSEMBLY_FAI" >> "${EXCLUSION_DIR}/foreground_and_gaps.bed"
sort -k1,1V -k2,2n "${EXCLUSION_DIR}/foreground_and_gaps.bed" | \
    bedtools merge -i - > "$EXCLUSION_MASK"

log "Exclusion mask: $(wc -l < $EXCLUSION_MASK) regions, total excluded: $(awk '{sum+=$3-$2}END{print sum}' $EXCLUSION_MASK) bp"

# ============================================================================
# SECTION 9: Mappability estimation
# ============================================================================
log "=== SECTION 9: Mappability estimation ==="

# Use H3K27ac (most broadly distributed mark) fragment coverage as mappability proxy
H3K27AC_FRAG="${FRAGMENT_DIR}/XG_152_fragments.bed"
MAPPABILITY_BG="${MAPPABILITY_DIR}/h3k27ac_coverage.bg"

if [[ -f "$MAPPABILITY_BG" ]] && [[ -s "$MAPPABILITY_BG" ]]; then
    log "SKIP mappability BedGraph: exists"
else
    log "Computing H3K27ac fragment coverage (mappability proxy)..."
    bedtools genomecov -i "$H3K27AC_FRAG" -g "${DATA_DIR}/chrom_sizes.txt" -bg > "$MAPPABILITY_BG"
    log "Mappability BedGraph: $(wc -l < $MAPPABILITY_BG) intervals"
fi

# Compute per-interval mappability scores for Set A
MAPPABILITY_SCORES="${MAPPABILITY_DIR}/setA_mappability_scores.txt"
if [[ -f "$MAPPABILITY_SCORES" ]] && [[ -s "$MAPPABILITY_SCORES" ]]; then
    log "SKIP mappability scores: exists"
else
    log "Computing per-interval mappability scores..."
    # Count fragments per interval, then compute mappable fraction
    # Mappable fraction = fraction of bases with >= 1 H3K27ac fragment
    bedtools coverage -a "${FOREGROUND_DIR}/setA_centroAnno_strict.bed" \
        -b "$H3K27AC_FRAG" > "${MAPPABILITY_DIR}/setA_coverage_raw.txt"

    awk -v OFS='\t' '{
        len = $3 - $2
        if (len > 0) {
            mappable = ($7 > 0 ? $7 : 0)
            fraction = mappable / len
            print $1, $2, $3, $4, fraction
        }
    }' "${MAPPABILITY_DIR}/setA_coverage_raw.txt" > "$MAPPABILITY_SCORES"

    log "Mappability scores: $(wc -l < $MAPPABILITY_SCORES) intervals"
    log "  Mean mappability: $(awk '{sum+=$5; n++} END{printf "%.4f", sum/n}' $MAPPABILITY_SCORES)"
    log "  Median mappability: $(awk '{a[NR]=$5} END{n=asort(a); print a[int(n/2)]}' $MAPPABILITY_SCORES)"
fi

# ============================================================================
# SECTION 10: QC summary
# ============================================================================
log "=== SECTION 10: QC summary ==="

cat > "${QC_DIR}/validation_summary.txt" << EOF
CENPA CUT&Tag Enrichment — Validation Summary
==============================================
Date: $(date)
Assembly: $ASSEMBLY_FASTA

Samples:
EOF

for sample in "${BAM_LIST[@]}"; do
    n_reads=$(grep "^$sample " "${QC_DIR}/library_sizes_reads.txt" | awk '{print $2}')
    n_frags=$(grep "^$sample " "${QC_DIR}/library_sizes_fragments.txt" | awk '{print $2}')
    echo "  $sample (${SAMPLE_ANTIBODY[$sample]}, ${SAMPLE_REPLICATE[$sample]}): $n_reads reads, $n_frags fragments" >> "${QC_DIR}/validation_summary.txt"
done

cat >> "${QC_DIR}/validation_summary.txt" << EOF

Foreground sets:
  Set A (strict): $(wc -l < ${FOREGROUND_DIR}/setA_centroAnno_strict.bed) intervals
  Set B (d=${MERGE_DIST_1}/pad=${FLANK_PAD_1}): $(wc -l < ${FOREGROUND_DIR}/setB_merged_d${MERGE_DIST_1}_pad${FLANK_PAD_1}.bed) regions
  Set B (d=${MERGE_DIST_2}/pad=${FLANK_PAD_2}): $(wc -l < ${FOREGROUND_DIR}/setB_merged_d${MERGE_DIST_2}_pad${FLANK_PAD_2}.bed) regions
  Set B (d=${MERGE_DIST_3}/pad=${FLANK_PAD_3}): $(wc -l < ${FOREGROUND_DIR}/setB_merged_d${MERGE_DIST_3}_pad${FLANK_PAD_3}.bed) regions
  Set C (HiCAT chr4): ${n_hicat_regions} occupied regions

Exclusion mask: $(wc -l < $EXCLUSION_MASK) regions
Gap regions: $(wc -l < $GAP_BED)
Chromosomes: $(wc -l < ${DATA_DIR}/chrom_sizes.txt)

Per-chromosome centroAnno intervals (chromosomes with <30):
$(awk '$1<30' ${QC_DIR}/setA_per_chromosome_counts.txt)

Mappability (H3K27ac proxy):
  Mean: $(awk '{sum+=$5; n++} END{printf "%.4f", sum/n}' $MAPPABILITY_SCORES)
  Median: $(awk '{a[NR]=$5} END{n=asort(a); print a[int(n/2)]}' $MAPPABILITY_SCORES)
EOF

log "QC summary written to ${QC_DIR}/validation_summary.txt"
log "=== 1_validate_and_prepare.sh DONE ==="
