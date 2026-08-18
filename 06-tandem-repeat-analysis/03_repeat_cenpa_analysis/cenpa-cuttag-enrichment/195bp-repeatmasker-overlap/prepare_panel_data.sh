#!/bin/bash
# ============================================================================
# prepare_panel_data.sh — build the 4 cached data tables for the
#   "195 bp repeats are L1" 3-panel supplementary figure.
#
#   panel_a_composition.tsv  — RepeatMasker class composition of the bin4
#                              (195 bp) vs bin6 (349 bp) arrays, as an
#                              EXCLUSIVE partition (each base counted once;
#                              L1 priority, then Unknown/Other on the non-L1
#                              remainder, remainder = Unannotated).
#   panel_b_consensus.tsv    — L1 consensus 5' 195 bp tandem-repeat blocks
#                              (from TRF on the consensus, already verified).
#   panel_c_locus.tsv        — RepeatMasker annotations across the example
#                              locus chr21:83,255,000-83,340,000 (the longest
#                              195 bp array), colored by family/class.
#   panel_c_arrays.tsv       — bin4 (195 bp) arrays in that window.
#
# Plot-only script reads these; re-run this only when data changes.
# ============================================================================
set -euo pipefail

BT=/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/genome-assembly-assessment/bin/bedtools
WD="$(cd "$(dirname "$0")" && pwd)"
SC=/tscc/lustre/ddn/scratch/jhc103
FAM=$SC/rm_hifiasm041425.family.bed
RES=$WD/results
mkdir -p "$RES" "$SC/rmcats"
cd "$WD"

# ── 1. Per-category RepeatMasker element files (from family bed) ────────────
awk -F'\t' '$5=="LINE/L1" && $4=="rnd-1_family-189"{print $1"\t"$2"\t"$3}' "$FAM" | sort -k1,1 -k2,2n > "$SC/rmcats/l1_fam189.bed"
awk -F'\t' '$5=="LINE/L1" && $4=="rnd-1_family-18"{print  $1"\t"$2"\t"$3}' "$FAM" | sort -k1,1 -k2,2n > "$SC/rmcats/l1_fam18.bed"
awk -F'\t' '$5=="LINE/L1" && $4!="rnd-1_family-189" && $4!="rnd-1_family-18"{print $1"\t"$2"\t"$3}' "$FAM" | sort -k1,1 -k2,2n > "$SC/rmcats/l1_other.bed"
awk -F'\t' '$5=="Unknown"{print $1"\t"$2"\t"$3}' "$FAM" | sort -k1,1 -k2,2n > "$SC/rmcats/unknown.bed"
awk -F'\t' '$5!="LINE/L1" && $5!="Unknown"{print $1"\t"$2"\t"$3}' "$FAM" | sort -k1,1 -k2,2n > "$SC/rmcats/other.bed"
cat "$SC/rmcats/l1_fam189.bed" "$SC/rmcats/l1_fam18.bed" "$SC/rmcats/l1_other.bed" \
  | sort -k1,1 -k2,2n | "$BT" merge -i - -d 0 > "$SC/rmcats/l1_union.bed"

# ── 2. Helper: union coverage of a category file clipped to target intervals ─
# cov_in <category.bed> <target.bed>  → echoes bp of category∩target (each base once)
cov_in () {
  "$BT" intersect -a "$2" -b "$1" -wao \
  | awk -F'\t' '$NF>0 { s=($2>$6?$2:$6); e=($3<$7?$3:$7); print $1"\t"s"\t"e }' \
  | sort -k1,1 -k2,2n | "$BT" merge -i - -d 0 | awk '{s+=$3-$2} END{print s+0}'
}

# ── 3. Panel (a): exclusive composition table ────────────────────────────────
partition () {  # $1 = arrays.bed ; writes "<array_group>\t<category>\t<bp>\t<pct>"
    local ARR="$1" TOTAL l1_ov
    TOTAL=$(awk '{s+=$3-$2} END{print s}' "$ARR")
    # L1 union coverage clipped to arrays → the L1-covered part
    l1=$(cov_in "$SC/rmcats/l1_union.bed" "$ARR")
    # exact L1-covered intervals (for subtract)
    "$BT" intersect -a "$ARR" -b "$SC/rmcats/l1_union.bed" -wao \
      | awk -F'\t' '$NF>0 { s=($2>$6?$2:$6); e=($3<$7?$3:$7); print $1"\t"s"\t"e }' \
      | sort -k1,1 -k2,2n | "$BT" merge -i - -d 0 > "$SC/rmcats/l1_clip.bed"
    # non-L1 remainder of the arrays
    "$BT" subtract -a "$ARR" -b "$SC/rmcats/l1_clip.bed" > "$SC/rmcats/nonl1.bed"
    NONL1=$(awk '{s+=$3-$2} END{print s+0}' "$SC/rmcats/nonl1.bed")
    UNK=$(cov_in "$SC/rmcats/unknown.bed" "$SC/rmcats/nonl1.bed")
    OTH=$(cov_in "$SC/rmcats/other.bed"  "$SC/rmcats/nonl1.bed")
    F189=$(cov_in "$SC/rmcats/l1_fam189.bed" "$ARR")
    F18=$(cov_in "$SC/rmcats/l1_fam18.bed"  "$ARR")
    FOTH=$(cov_in "$SC/rmcats/l1_other.bed" "$ARR")
    UNANN=$(( NONL1 - UNK - OTH ))
    # (l1 shown as the three family pieces; unannotated absorbs any tiny gap)
    printf "%s\tl1_fam189\t%d\t%.4f\n"  "$LBL" "$F189" "$(awk "BEGIN{print 100*$F189/$TOTAL}")"
    printf "%s\tl1_fam18\t%d\t%.4f\n"   "$LBL" "$F18"  "$(awk "BEGIN{print 100*$F18/$TOTAL}")"
    printf "%s\tl1_other\t%d\t%.4f\n"   "$LBL" "$FOTH" "$(awk "BEGIN{print 100*$FOTH/$TOTAL}")"
    printf "%s\tunknown\t%d\t%.4f\n"    "$LBL" "$UNK"  "$(awk "BEGIN{print 100*$UNK/$TOTAL}")"
    printf "%s\tother\t%d\t%.4f\n"      "$LBL" "$OTH"  "$(awk "BEGIN{print 100*$OTH/$TOTAL}")"
    printf "%s\tunannotated\t%d\t%.4f\n" "$LBL" "$UNANN" "$(awk "BEGIN{print 100*$UNANN/$TOTAL}")"
}

{
  echo -e "array_group\tcategory\tbp\tpct"
  LBL=bin4_195bp; partition bin4_arrays_seq.bed
  LBL=bin6_349bp; partition bin6_arrays_seq.bed
} > "$RES/panel_a_composition.tsv"
echo "panel_a written:"; cat "$RES/panel_a_composition.tsv"

# ── 4. Panel (b): L1 consensus 5' repeat blocks (verified by TRF) ──────────
{
  echo -e "family\tconsensus_len\trepeat_start\trepeat_end\tunit_len\tblocks_csv"
  echo -e "rnd-1_family-189\t7758\t1\t804\t195\t1-195|196-390|391-585|586-780|781-804"
  echo -e "rnd-1_family-18\t4933\t1\t438\t195\t1-195|196-390|391-438"
} > "$RES/panel_b_consensus.tsv"
echo "panel_b written"

# ── 5. Panel (c): example locus chr21 (seq22):83,255,000-83,340,000 ────────
W0=83255000; W1=83340000
awk -F'\t' -v a=$W0 -v b=$W1 '$1=="seq22" && $2 < b && $3 > a {
  cat = "Other";
  if ($5 == "LINE/L1")      cat = ($4=="rnd-1_family-189" ? "L1 fam-189" : ($4=="rnd-1_family-18" ? "L1 fam-18" : "L1 other"));
  else if ($5 == "Unknown") cat = "Unknown";
  printf "chr21\t%d\t%d\t%s\t%s\t%s\n", $2, $3, cat, $4, $6
}' "$FAM" | sort -k2,2n > "$RES/panel_c_locus.tsv"
awk -F'\t' -v a=$W0 -v b=$W1 '$1=="seq22" && $2 < b && $3 > a {
  printf "chr21\t%d\t%d\t%s\n", $2, $3, $4
}' bin4_arrays_seq.bed | sort -k2,2n > "$RES/panel_c_arrays.tsv"
echo "panel_c written: $(wc -l < "$RES/panel_c_locus.tsv") RM annotations, $(wc -l < "$RES/panel_c_arrays.tsv") bin4 arrays"
head -25 "$RES/panel_c_locus.tsv"
echo "--- bin4 arrays in window:"
cat "$RES/panel_c_arrays.tsv"
