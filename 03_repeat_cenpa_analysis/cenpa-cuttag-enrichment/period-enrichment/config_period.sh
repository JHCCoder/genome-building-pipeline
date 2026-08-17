#!/bin/bash
# config_period.sh — paths and parameters for TRF period-size stratified enrichment
# Sources the parent config.sh for shared paths

set -euo

PARENT_CONFIG="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/config.sh"
if [[ -f "$PARENT_CONFIG" ]]; then
    source "$PARENT_CONFIG"
else
    echo "ERROR: Parent config not found: $PARENT_CONFIG"
    exit 1
fi

# ============================================================================
# Period-enrichment specific paths
# ============================================================================
PERIOD_DIR="${FIGURE_DIR}/period-enrichment"
PERIOD_DATA_DIR="${PERIOD_DIR}/data"
PERIOD_SCRIPTS_DIR="${PERIOD_DIR}/scripts"
PERIOD_PLOTS_DIR="${PERIOD_DIR}/plots"
PERIOD_RESULTS_DIR="${PERIOD_DIR}/results"

# TRF input
TRF_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/genome-annotation/trf-tandem-repeat"
TRF_REPEAT_DF="${TRF_DIR}/repeat_df_degu.tsv"

# Output files
TRF_BED_ALL="${PERIOD_DATA_DIR}/trf_all_period_bins.bed"
TRF_BED_CHR="${PERIOD_DATA_DIR}/trf_chr_only_period_bins.bed"
BIN_STATS="${PERIOD_DATA_DIR}/period_bin_statistics.csv"

# Signal extraction output
COUNTS_DIR="${PERIOD_DATA_DIR}/counts"

# Coverage files
COVERAGE_FILES=(
    "${COVERAGE_DIR}/XG_150_perbase.txt.gz"
    "${COVERAGE_DIR}/XG_151_perbase.txt.gz"
    "${COVERAGE_DIR}/XG_152_perbase.txt.gz"
    "${COVERAGE_DIR}/XG_153_perbase.txt.gz"
)

# Samples
SAMPLES=("XG_150" "XG_151" "XG_152" "XG_153")

# Chromosomes (chr1-28, chrX, chrY)
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chrX chrY"

# ============================================================================
# Slurm settings
# ============================================================================
SLURM_PARTITION="condo"
SLURM_QOS="condo"
SLURM_ACCOUNT="csd788"
CLUSTER_LOG_DIR="/tscc/nfs/home/jhc103/cluster-logs"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

init_period_dirs() {
    mkdir -p "$PERIOD_DATA_DIR" "$PERIOD_SCRIPTS_DIR" "$PERIOD_PLOTS_DIR" \
             "$PERIOD_RESULTS_DIR" "$COUNTS_DIR"
}

log "config_period.sh loaded successfully"
