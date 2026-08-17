#!/bin/bash
# config_349bp.sh — paths and parameters for 349-bp repeat CENP-A enrichment analysis
# Sources the parent config.sh for shared paths (BAMs, assembly, fragment BEDs, etc.)

set -euo

# Source parent config
PARENT_CONFIG="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/config.sh"
if [[ -f "$PARENT_CONFIG" ]]; then
    source "$PARENT_CONFIG"
else
    echo "ERROR: Parent config not found: $PARENT_CONFIG"
    exit 1
fi

# ============================================================================
# 349-bp specific paths
# ============================================================================
ANALYSIS_349_DIR="${FIGURE_DIR}/349-bp"
DATA_349_DIR="${ANALYSIS_349_DIR}/data"
FOREGROUND_349_DIR="${DATA_349_DIR}/foregrounds"
BACKGROUND_349_DIR="${DATA_349_DIR}/backgrounds"
COUNTS_349_DIR="${DATA_349_DIR}/counts"
PLOTS_349_DIR="${ANALYSIS_349_DIR}/plots"
RESULTS_349_DIR="${ANALYSIS_349_DIR}/results"

# ============================================================================
# TRF input data
# ============================================================================
TRF_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/code/command-line-script/genome-annotation/trf-tandem-repeat"
TRF_REPEAT_DF="${TRF_DIR}/repeat_df_degu.tsv"

# ============================================================================
# 349-bp parameters
# ============================================================================
PERIOD_MIN=345
PERIOD_MAX=355
MERGE_DIST=50000          # merge nearby 349-bp repeats within 50 kb
N_SHUFFLES=1000           # background iterations
SHUFFLE_SEED=20260722     # different seed from V2 for independence

# Chromosomes: chr1-28, chrX only (chrY has no 349-bp TRF repeats)
CHROMOSOMES_349="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chrX"

# ============================================================================
# Slurm settings
# ============================================================================
SLURM_PARTITION="condo"
SLURM_QOS="condo"
SLURM_ACCOUNT="csd788"
CLUSTER_LOG_DIR="/tscc/nfs/home/jhc103/cluster-logs"

# ============================================================================
# Output file paths
# ============================================================================
TRF_349BP_RAW="${FOREGROUND_349_DIR}/trf_349bp_raw.bed"
HORS_349BP_RAW="${FOREGROUND_349_DIR}/hors_349bp_raw.bed"
FOREGROUND_MERGED="${FOREGROUND_349_DIR}/349bp_merged.bed"
BG_SHUFFLED="${BACKGROUND_349_DIR}/bg_shuffled.bed"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

init_349_dirs() {
    mkdir -p "$FOREGROUND_349_DIR" "$BACKGROUND_349_DIR" "$COUNTS_349_DIR" \
             "$PLOTS_349_DIR" "$RESULTS_349_DIR"
}

log "config_349bp.sh loaded successfully"
