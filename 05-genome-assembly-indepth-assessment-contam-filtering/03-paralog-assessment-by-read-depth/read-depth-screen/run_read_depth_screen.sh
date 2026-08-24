#!/bin/bash
#SBATCH -J read_depth_screen
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH -t 2-00:00:00
#SBATCH --mem=256G
# --- Cluster-specific (Slurm on TSCC/UCSD) — adjust for your scheduler ---
#SBATCH -o /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.out
#SBATCH -e /tscc/nfs/home/jhc103/cluster-logs/%x.%j.%N.err
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd788
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user=you@example.com   # TODO: your email

# This script needs TWO directories, and they are not the same:
#
#   WORK_DIR  — data: the HiFi BAM, ../paralog_families.tsv, and all outputs.
#               The working copy of this analysis lives under figure/, outside
#               the repo.
#   CODE_DIR_SCREEN — code: read_depth_screen.py, which sits next to this script
#               inside the repo.
#
# Under Slurm, BASH_SOURCE[0] points at the *spool copy* of the script
# (/cm/local/.../slurm_script), so it can supply neither. WORK_DIR comes from
# SLURM_SUBMIT_DIR (where sbatch was invoked); CODE_DIR_SCREEN is derived from
# the location of the config.sh we load. Both fall back sensibly for a plain
# `bash run_read_depth_screen.sh`.
if [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    WORK_DIR="$SLURM_SUBMIT_DIR"
else
    WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
fi

# Load shared configuration (config.sh at the repo root).
# The run directory is not necessarily inside the repo, so honour an explicit
# CONFIG_SH first and only then walk upwards.
if [[ -n "${CONFIG_SH:-}" && -f "$CONFIG_SH" ]]; then
    _config_sh="$CONFIG_SH"
else
    _repo_root="$WORK_DIR"
    while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
    if [[ ! -f "$_repo_root/config.sh" ]]; then
        echo "ERROR: could not locate config.sh above $WORK_DIR." >&2
        echo "       Set CONFIG_SH=/path/to/config.sh (sbatch --export=ALL,CONFIG_SH=...)." >&2
        exit 1
    fi
    _config_sh="$_repo_root/config.sh"
    unset _repo_root
fi
source "$_config_sh"

CODE_DIR_SCREEN="$(cd "$(dirname "$_config_sh")" && pwd)/05-genome-assembly-indepth-assessment-contam-filtering/03-paralog-assessment-by-read-depth/read-depth-screen"
unset _config_sh
if [[ ! -f "$CODE_DIR_SCREEN/read_depth_screen.py" ]]; then
    echo "ERROR: read_depth_screen.py not found in $CODE_DIR_SCREEN" >&2
    exit 1
fi

source ~/.bashrc
# Temporarily disable -u so conda's openjdk_activate.sh doesn't crash
set +u
conda activate "$ENV_ASSESSMENT"
set -euo pipefail

ASSEMBLY="$READ_DEPTH_ASSEMBLY"
HIFI_READS="$READ_DEPTH_HIFI_DIR"
PARALOG_FAMILIES="$WORK_DIR/../paralog_families.tsv"   # from build_paralog_families.py
REPEAT_BED="$PARALOG_REPEAT_BED"
GENOME_LENGTH_FILE="$PARALOG_GENOME_LENGTHS"

BAM_PREFIX="hifi_male_map-hifi_contamFiltered_chrAssigned"
OUT_DIR="${WORK_DIR}"
mkdir -p "${OUT_DIR}"

# Step 1: map HiFi reads with minimap2 map-hifi → sort → index
echo "STEP 1: mapping HiFi reads with minimap2 map-hifi"
BAM_FILE="${OUT_DIR}/${BAM_PREFIX}.bam"
if [ -f "${BAM_FILE}" ] && [ -f "${BAM_FILE}.bai" ]; then
    echo "BAM and index already exist: ${BAM_FILE} — skipping mapping."
else
    READ_FILES=()
    for f in "${HIFI_READS}"/*.fastq.gz; do
        [ -f "$f" ] && READ_FILES+=("$f")
    done
    echo "Found ${#READ_FILES[@]} HiFi read files."
    # -F 0x900 excludes supplementary (0x800) and secondary (0x100) alignments
    minimap2 -ax map-hifi -t 36 "${ASSEMBLY}" "${READ_FILES[@]}" \
        | samtools view -@ 8 -F 0x900 -Sb \
        | samtools sort -@ 8 -o "${BAM_FILE}" \
        && samtools index -@ 8 "${BAM_FILE}"
fi

# Step 2: coverage analysis / haplotig screening
#
# read_depth_screen.py caches each phase by checking whether its output file
# already exists. That is a real hazard after the gene set changes: the run
# loads the NEW paralog_families.tsv but then classifies using a cached
# per_gene_coverage.tsv / pair_classification.tsv built from the OLD gene set,
# finishes in under a minute, and reports success with stale results.
#
# So pass READ_DEPTH_EXTRA_ARGS="--no-cache" whenever the annotation or the
# paralog families have changed:
#
#   READ_DEPTH_EXTRA_ARGS="--no-cache" sbatch --export=ALL,... run_read_depth_screen.sh
#
# Left empty by default so routine re-runs still reuse the expensive
# genome-wide coverage pass.
echo "STEP 2: coverage analysis and haplotig screening"
if [[ -n "${READ_DEPTH_EXTRA_ARGS:-}" ]]; then
    echo "  extra args: ${READ_DEPTH_EXTRA_ARGS}"
fi
python3 "${CODE_DIR_SCREEN}/read_depth_screen.py" \
    --bam "${BAM_FILE}" \
    --fasta "${ASSEMBLY}" \
    --paralog-families "${PARALOG_FAMILIES}" \
    --genome-lengths "${GENOME_LENGTH_FILE}" \
    --out-dir "${OUT_DIR}" \
    --repeat-bed "${REPEAT_BED}" \
    --threads 24 \
    --flank-bp 25000 \
    ${READ_DEPTH_EXTRA_ARGS:-}

echo "PIPELINE COMPLETE — output: ${OUT_DIR}"
