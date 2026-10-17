#!/bin/bash
# Step 1: collect the peak–gene linked peaks from the four OctDeg2.0 Cell Ranger ARC reruns
# (Liftoff-only reference) and split them by whether their OctDeg2.0 sequence is covered by
# the OctDeg1.0 <-> OctDeg2.0 minimap2 alignment used in Figure 3A ($FIG3_ALIGNED_BED).
#
#   shared candidate : >= 50% of the peak lies in aligned OctDeg2.0 sequence
#   unaligned        : <  50% aligned -> candidate "new" cCRE, checked against OctDeg1.0
#                      sequence directly in step 2 (minimap2 misses collapsed and high-copy
#                      repeats, so "unaligned" overstates "absent from OctDeg1.0").
#
# Only rows of type "peak-gene" are used (each peak–gene pair also appears as "gene-peak";
# peak–peak co-accessibility links are excluded).
# Lightweight (bedtools on ~10^5 intervals, < 1 min); runs on a login node.
#
# Outputs: $FIG3D_CCRE_DIR/01_*  and  results/01_alignment_overlap_summary.tsv
set -euo pipefail

_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
_repo_root="$_script_dir"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

BEDTOOLS="$BEDTOOLS_BIN"
SAMPLES=(181_PFC 181_dHIP 6997_dHIP 7000_PFC)
OUT="$FIG3D_CCRE_DIR"
RES="$_script_dir/results"
mkdir -p "$OUT" "$RES"

# Merged OctDeg2.0 intervals covered by any OctDeg1.0 alignment
sort -k1,1 -k2,2n "$FIG3_ALIGNED_BED" | "$BEDTOOLS" merge -i - > "$OUT/01_octdeg2_aligned_merged.bed"

summary="$RES/01_alignment_overlap_summary.tsv"
echo -e "set\tlinked_peaks\tshared_candidates\tunaligned_candidates\tpct_unaligned\tfully_unaligned" > "$summary"

classify () {   # $1 = label, $2 = bed (chr start end ...)
  "$BEDTOOLS" coverage -a "$2" -b "$OUT/01_octdeg2_aligned_merged.bed" > "$OUT/01_${1}.alignedFrac.bed"
  awk -F'\t' -v s="$1" '{f=$NF; n++; if (f>=0.5) sh++; else un++; if (f==0) z++}
       END {printf "%s\t%d\t%d\t%d\t%.2f\t%d\n", s, n, sh, un, 100*un/n, z}' \
       "$OUT/01_${1}.alignedFrac.bed" >> "$summary"
}

# Per-sample unique linked peaks (chr, start, end, sample, linked genes)
for s in "${SAMPLES[@]}"; do
  bedpe="$FIG3D_CR_DIR/count/${s}_OctDeg2_liftoffOnly/outs/analysis/feature_linkage/feature_linkage.bedpe"
  awk -F'\t' -v s="$s" 'BEGIN{OFS="\t"} $13=="peak-gene" {k=$1"\t"$2"\t"$3; g[k]=(k in g ? g[k]"," : "") $7}
       END {for (k in g) print k, s, g[k]}' "$bedpe" | sort -k1,1 -k2,2n > "$OUT/01_linked_peaks_${s}.bed"
  classify "$s" "$OUT/01_linked_peaks_${s}.bed"
done

# Union across samples: overlapping peaks merged, with the samples that called them
cat "$OUT"/01_linked_peaks_{181_PFC,181_dHIP,6997_dHIP,7000_PFC}.bed | sort -k1,1 -k2,2n \
  | "$BEDTOOLS" merge -i - -c 4,4 -o count_distinct,distinct > "$OUT/01_linked_peaks_union.bed"
classify union "$OUT/01_linked_peaks_union.bed"

# Candidate "new" peaks (< 50% aligned) for the direct OctDeg1.0 check in step 2
awk -F'\t' 'BEGIN{OFS="\t"} $NF<0.5 {print $1,$2,$3,$1":"$2"-"$3,$4,$5,$NF}' \
  "$OUT/01_union.alignedFrac.bed" > "$OUT/01_union_unaligned_candidates.bed"
"$BEDTOOLS" getfasta -fi "$MITO_ASSEMBLY" -bed "$OUT/01_union_unaligned_candidates.bed" -nameOnly \
  > "$OUT/01_union_unaligned_candidates.fa"

column -t "$summary"
echo "candidates: $(wc -l < "$OUT/01_union_unaligned_candidates.bed") peaks -> $OUT/01_union_unaligned_candidates.fa"
