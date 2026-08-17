#!/bin/bash
# 00_discover_inputs.sh — Discover and inspect input files for purge_dups vs BISER enrichment analysis
# Run: bash scripts/00_discover_inputs.sh

set -euo pipefail

PROJECT_ROOT="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
PURGE_ROOT="${PROJECT_ROOT}/output/outputs-from-purge-duplicate/hifiasm-041425-scaffolded-assembly"
BISER_ROOT="${PROJECT_ROOT}/code/command-line-script/genome-annotation/biser/hifiasm-041425"
OUTPUT_ROOT="${PROJECT_ROOT}/figure/segdup-purgeDup-overlap/purge-dups-biser-enrichment"
LOG_DIR="${OUTPUT_ROOT}/logs"

mkdir -p "${LOG_DIR}"

REPORT="${LOG_DIR}/input_discovery_report.txt"

{
echo "============================================"
echo "INPUT DISCOVERY REPORT"
echo "Generated: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "============================================"
echo ""
echo "PROJECT_ROOT: ${PROJECT_ROOT}"
echo "PURGE_ROOT:   ${PURGE_ROOT}"
echo "BISER_ROOT:   ${BISER_ROOT}"
echo ""

# ---- Section 1: purge_dups directory ----
echo "============================================"
echo "SECTION 1: purge_dups DIRECTORY CONTENTS"
echo "============================================"
echo ""
echo "--- All files ---"
find "$PURGE_ROOT" -maxdepth 6 -type f | sort
echo ""
echo "--- File sizes (bytes) ---"
find "$PURGE_ROOT" -maxdepth 6 -type f -printf '%s\t%p\n' | sort -n
echo ""

echo "--- Searching for category labels ---"
grep -RIlE 'HAPLOTIG|HIGHCOV|JUNK|OVLP|REPEAT' "$PURGE_ROOT" 2>/dev/null || echo "(none found)"
echo ""

echo "--- dups.bed: first 30 lines ---"
head -30 "${PURGE_ROOT}/dups.bed" 2>/dev/null || echo "(file not found)"
echo ""

echo "--- dups.bed: column counts ---"
awk '{print NF}' "${PURGE_ROOT}/dups.bed" 2>/dev/null | sort | uniq -c || echo "(file not found)"
echo ""

echo "--- dups.bed: category counts ---"
grep -oE 'HAPLOTIG|HIGHCOV|JUNK|OVLP|REPEAT' "${PURGE_ROOT}/dups.bed" 2>/dev/null | sort | uniq -c || echo "(no categories found)"
echo ""

echo "--- dups.bed: total line count ---"
wc -l "${PURGE_ROOT}/dups.bed" 2>/dev/null || echo "(file not found)"
echo ""

echo "--- dups.3col.bed: first 10 lines ---"
head -10 "${PURGE_ROOT}/dups.3col.bed" 2>/dev/null || echo "(file not found)"
echo ""

echo "--- dups_haplo.bed: full contents ---"
cat "${PURGE_ROOT}/dups_haplo.bed" 2>/dev/null || echo "(file not found)"
echo ""

echo "--- PB.stat ---"
cat "${PURGE_ROOT}/PB.stat" 2>/dev/null || echo "(file not found)"
echo ""

echo "--- calcults.log ---"
cat "${PURGE_ROOT}/calcults.log" 2>/dev/null || echo "(file not found)"
echo ""

echo "--- cutoffs ---"
cat "${PURGE_ROOT}/cutoffs" 2>/dev/null || echo "(file not found)"
echo ""

# ---- Section 2: BISER directory ----
echo ""
echo "============================================"
echo "SECTION 2: BISER DIRECTORY CONTENTS"
echo "============================================"
echo ""
echo "--- All files ---"
find "$BISER_ROOT" -maxdepth 6 -type f | sort
echo ""
echo "--- File sizes (bytes) ---"
find "$BISER_ROOT" -maxdepth 6 -type f -printf '%s\t%p\n' | sort -n
echo ""

echo "--- segdup_output_mod.bedpe: first 20 lines ---"
head -20 "${BISER_ROOT}/segdup_output_mod.bedpe" 2>/dev/null || echo "(file not found)"
echo ""
echo "--- segdup_output_mod.bedpe: line count ---"
wc -l "${BISER_ROOT}/segdup_output_mod.bedpe" 2>/dev/null || echo "(file not found)"
echo ""
echo "--- segdup_output_mod.bedpe: column count ---"
head -5 "${BISER_ROOT}/segdup_output_mod.bedpe" 2>/dev/null | awk '{print NF}' | sort | uniq -c || echo "(file not found)"

echo ""
echo "--- segdup_output.bedpe: first 5 lines ---"
head -5 "${BISER_ROOT}/segdup_output.bedpe" 2>/dev/null || echo "(file not found)"
echo ""
echo "--- segdup_output.bedpe: line count ---"
wc -l "${BISER_ROOT}/segdup_output.bedpe" 2>/dev/null || echo "(file not found)"
echo ""
echo "--- segdup_output.bedpe: column count ---"
head -5 "${BISER_ROOT}/segdup_output.bedpe" 2>/dev/null | awk '{print NF}' | sort | uniq -c || echo "(file not found)"

echo ""
echo "--- segdup_output_duplicateLinkRemoved.bedpe: first 5 lines ---"
head -5 "${BISER_ROOT}/segdup_output_duplicateLinkRemoved.bedpe" 2>/dev/null || echo "(file not found)"
echo ""
echo "--- segdup_output_duplicateLinkRemoved.bedpe: line count ---"
wc -l "${BISER_ROOT}/segdup_output_duplicateLinkRemoved.bedpe" 2>/dev/null || echo "(file not found)"

echo ""
echo "--- segdup_output_hifiasm-041425.elem.txt: first 10 lines ---"
head -10 "${BISER_ROOT}/segdup_output_hifiasm-041425.elem.txt" 2>/dev/null || echo "(file not found)"
echo ""
echo "--- segdup_output_hifiasm-041425.elem.txt: line count ---"
wc -l "${BISER_ROOT}/segdup_output_hifiasm-041425.elem.txt" 2>/dev/null || echo "(file not found)"

echo ""
echo "--- Chromosomes in segdup_output_mod.bedpe (arm1) ---"
cut -f1 "${BISER_ROOT}/segdup_output_mod.bedpe" 2>/dev/null | sort -V | uniq -c | sort -k2 -V
echo ""
echo "--- Chromosomes in segdup_output_mod.bedpe (both arms, sorted) ---"
cat <(cut -f1 "${BISER_ROOT}/segdup_output_mod.bedpe") <(cut -f4 "${BISER_ROOT}/segdup_output_mod.bedpe") 2>/dev/null | sort -V | uniq -c | sort -k2 -V
echo ""
echo "--- Unique chromosome count ---"
cat <(cut -f1 "${BISER_ROOT}/segdup_output_mod.bedpe") <(cut -f4 "${BISER_ROOT}/segdup_output_mod.bedpe") 2>/dev/null | sort -u | wc -l

echo ""
echo "--- Unplaced scaffolds (seq*) in BISER ---"
cat <(cut -f1 "${BISER_ROOT}/segdup_output_mod.bedpe") <(cut -f4 "${BISER_ROOT}/segdup_output_mod.bedpe") 2>/dev/null | sort -u | grep -E '^seq'

# ---- Section 3: Assembly FASTA ----
echo ""
echo "============================================"
echo "SECTION 3: ASSEMBLY FASTA"
echo "============================================"
echo ""

GENOME_BROWSER_FA="${PROJECT_ROOT}/degu-genome-browser-pythonVersion/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta"
echo "--- Genome browser assembly ---"
ls -la "${GENOME_BROWSER_FA}" 2>/dev/null || echo "(not found)"
echo ""
echo "--- Genome browser assembly .fai ---"
if [ -f "${GENOME_BROWSER_FA}.fai" ]; then
    cat "${GENOME_BROWSER_FA}.fai"
    echo ""
    echo "--- Total sequences in .fai ---"
    wc -l "${GENOME_BROWSER_FA}.fai"
    echo ""
    echo "--- Total assembly length ---"
    awk '{sum+=$2} END {print sum}' "${GENOME_BROWSER_FA}.fai"
else
    echo "(no .fai found)"
fi

echo ""
echo "--- BISER assembly FASTA ---"
BISER_ASM="${BISER_ROOT}/deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_041325.asm.hic.p_ctg.fa"
ls -la "${BISER_ASM}" 2>/dev/null || echo "(not found)"

# ---- Section 4: File selection recommendations ----
echo ""
echo "============================================"
echo "SECTION 4: FILE SELECTION"
echo "============================================"
echo ""
echo "SELECTED purge_dups file: ${PURGE_ROOT}/dups.bed"
echo "  Reason: Contains 5-column BED with category labels (HAPLOTIG, HIGHCOV, JUNK, OVLP, REPEAT)"
echo "  245 total intervals across all categories"
echo ""
echo "SELECTED BISER file: ${BISER_ROOT}/segdup_output_mod.bedpe"
echo "  Reason: Simplified 6-column BEDPE (chr1, start1, end1, chr2, start2, end2)"
echo "  208,156 duplication pairs, no CIGAR/metadata columns"
echo "  Alternative candidates:"
echo "    - segdup_output.bedpe: 14-column full format, same pairs"
echo "    - segdup_output_duplicateLinkRemoved.bedpe: 203,532 pairs, has extra index column"
echo ""
echo "SELECTED assembly: ${GENOME_BROWSER_FA}"
echo "  Reason: Has existing .fai, matches chromosome naming in both purge_dups and BISER"
echo "  Contains: chr1-chr28, chrX, chrY, seq31-seq288, NC_020661.1 (mito)"
echo ""
echo "BISER was run with --keep-contigs (or equivalent): seq scaffolds are present in output"

} > "${REPORT}"

echo "Discovery report written to: ${REPORT}"
echo ""
echo "=== Quick summary ==="
echo "purge_dups intervals: $(wc -l < "${PURGE_ROOT}/dups.bed")"
echo "BISER pairs: $(wc -l < "${BISER_ROOT}/segdup_output_mod.bedpe")"
echo "Assembly sequences: $(wc -l < "${GENOME_BROWSER_FA}.fai")"
