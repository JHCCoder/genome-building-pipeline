#!/bin/bash
# 10_call_cenpa_peaks.sh — call CENP-A peaks with MACS2 (mac2 env) from the
# CUT&Tag fragment BEDs, then redefine each chromosome's CENP-A core as the
# peak interval overlapping the strict-signal peak window.
#
# CUT&Tag fragment BEDs are 3-col `chr start end` (k=1 unique-mapped fragments).
# We use the standard CUT&Tag recipe: --nomodel --shift -100 --extsize 200,
# with H3K27ac (XG_152+153) as the negative control.

set -euo
MACS=/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/mac2/bin/macs2
WORK="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
SRC="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
OUT="${WORK}/data/macs2_peaks"
mkdir -p "$OUT"

FRAG150="${SRC}/data/fragments/XG_150_fragments.bed"
FRAG151="${SRC}/data/fragments/XG_151_fragments.bed"
CTRL152="${SRC}/data/fragments/XG_152_fragments.bed"
CTRL153="${SRC}/data/fragments/XG_153_fragments.bed"

echo "[$(date '+%H:%M:%S')] calling CENP-A peaks (CENP-A reps vs H3K27ac)..."
"$MACS" callpeak \
    -t "$FRAG150" "$FRAG151" \
    -c "$CTRL152" "$CTRL153" \
    -f BED -g 3381500000 \
    --nomodel --shift -100 --extsize 200 \
    --keep-dup all -q 0.05 \
    --outdir "$OUT" -n cenpa \
    --bdg --SPMR 2>&1 | tail -15

echo "[$(date '+%H:%M:%S')] peaks: $(wc -l < ${OUT}/cenpa_peaks.narrowPeak 2>/dev/null || echo 0)"
echo "DONE"
