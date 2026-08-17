#!/bin/bash
#SBATCH -J map_and_coverage
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 24:00:00
#SBATCH --mem=256G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END
#SBATCH --mail-user=you@example.com   # TODO: your email

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

# ====== RUN PARAMETERS (edit these for each run) ====== #
READ_TYPE="long"      # short (bwa-mem2 mem) | long (minimap2 map-hifi)
MULTIMAP="yes"        # yes = keep secondary alignments (multi-mapping) | no = unique (primary only)
FILTER="yes"          # yes = drop supplementary alignments + MAPQ filter | no = no filtering
MAPQ=30               # mapping-quality threshold, applied only when FILTER=yes
THREADS=16
WINDOW_SIZE=15000     # coverage window size (bp)
ASSEMBLY_ALIAS="hifiasm_041425_scaffolded_juiceBox_sorted_chrAssigned"   # label used in output names

# ====== build mode tag + flags ====== #
case "$READ_TYPE" in
  short|long) ;;
  *) echo "ERROR: READ_TYPE must be 'short' or 'long'"; exit 1 ;;
esac

mode_tag=""
if [ "$MULTIMAP" = "yes" ]; then mode_tag="multimap"; else mode_tag="unique"; fi
if [ "$FILTER"   = "yes" ]; then mode_tag="${mode_tag}_filtered"; else mode_tag="${mode_tag}_unfiltered"; fi

# SAM flags to drop (bitwise OR): secondary=0x100, supplementary=0x800
drop_flag=0
[ "$MULTIMAP" = "no"  ] && drop_flag=$((drop_flag | 0x100))
[ "$FILTER"   = "yes" ] && drop_flag=$((drop_flag | 0x800))

view_flags=(-@ "$THREADS")
[ "$drop_flag" -ne 0 ] && view_flags+=(-F "$drop_flag")
[ "$FILTER" = "yes" ]  && view_flags+=(-q "$MAPQ")

bam="${ASSEMBLY_ALIAS}_${READ_TYPE}_read_${mode_tag}_alignment_sorted.bam"

mkdir -p "$COVERAGE_OUT_DIR"
cd "$COVERAGE_OUT_DIR" || exit 1

# ====== mapping ====== #
source ~/.bashrc   # initialize conda (adjust for your setup)
if [ "$READ_TYPE" = "short" ]; then
    conda activate "$ENV_BULK_HIC"

    # bwa-mem2 index (created next to the assembly if absent)
    if [ ! -f "${COVERAGE_ASSEMBLY}.0123" ]; then
        bwa-mem2 index "$COVERAGE_ASSEMBLY" || exit 1
    fi

    bwa_flags=(-t "$THREADS")
    [ "$MULTIMAP" = "yes" ] && bwa_flags+=(-a)   # -a: report all alignments (multi-mapping)

    bwa-mem2 mem "${bwa_flags[@]}" "$COVERAGE_ASSEMBLY" \
        <(zcat $COVERAGE_SHORT_R1) \
        <(zcat $COVERAGE_SHORT_R2) | \
        samtools view "${view_flags[@]}" -Sb | \
        samtools sort -@ "$THREADS" -o "$bam" - || exit 1
else
    conda activate "$ENV_ASSESSMENT"

    mm2_flags=(-ax map-hifi -t "$THREADS")
    [ "$MULTIMAP" = "no" ] && mm2_flags+=(--secondary=no)   # suppress secondary alignments

    minimap2 "${mm2_flags[@]}" "$COVERAGE_ASSEMBLY" $COVERAGE_HIFI_READS | \
        samtools view "${view_flags[@]}" -Sb | \
        samtools sort -@ "$THREADS" -o "$bam" - || exit 1
fi
samtools index -@ "$THREADS" "$bam"

# ====== genome windows (assembly-derived, shared across modes) ====== #
if [ ! -s "${COVERAGE_ASSEMBLY}.fai" ]; then
    samtools faidx "$COVERAGE_ASSEMBLY"
fi
cut -f1,2 "${COVERAGE_ASSEMBLY}.fai" > "${ASSEMBLY_ALIAS}_genome_length.txt"

# ====== coverage (bedtools) ====== #
conda activate "$ENV_BEDTOOLS"

bedtools makewindows -g "${ASSEMBLY_ALIAS}_genome_length.txt" -w "$WINDOW_SIZE" > "${ASSEMBLY_ALIAS}_${WINDOW_SIZE}bp_windows.bed"
bedtools coverage -a "${ASSEMBLY_ALIAS}_${WINDOW_SIZE}bp_windows.bed" -b "$bam" -sorted -g "${ASSEMBLY_ALIAS}_genome_length.txt" -mean \
    > "coverage_${READ_TYPE}_read_${ASSEMBLY_ALIAS}_${mode_tag}_${WINDOW_SIZE}bp_windows.tsv"

cat "${ASSEMBLY_ALIAS}_genome_length.txt" | awk '{print $1 "\t0\t" $2}' > "${ASSEMBLY_ALIAS}_contigs.bed"
bedtools coverage -a "${ASSEMBLY_ALIAS}_contigs.bed" -b "$bam" -sorted -g "${ASSEMBLY_ALIAS}_genome_length.txt" -mean \
    > "coverage_${READ_TYPE}_read_${ASSEMBLY_ALIAS}_${mode_tag}_contig.tsv"

awk '{total_len += $3; total_cov += $3 * $4} END {print "Whole genome average coverage:", total_cov / total_len}' \
    "coverage_${READ_TYPE}_read_${ASSEMBLY_ALIAS}_${mode_tag}_contig.tsv" \
    > "coverage_${READ_TYPE}_read_${ASSEMBLY_ALIAS}_${mode_tag}_whole_assembly.tsv"

echo "SUCCESS: ${READ_TYPE} reads mapped (${mode_tag}) and coverage written to ${COVERAGE_OUT_DIR}"
