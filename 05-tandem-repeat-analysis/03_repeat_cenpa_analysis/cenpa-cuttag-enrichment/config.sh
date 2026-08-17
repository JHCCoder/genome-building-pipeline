#!/bin/bash
# config.sh — shared paths, parameters, and environment for CENPA enrichment analysis
# Source this file in all downstream scripts: source config.sh

set -euo
# NOTE: no pipefail — SIGPIPE from `| head`/`| tail` on samtools streams is
# expected behavior, not an error. With pipefail, every pipeline ending in head
# would kill the script.

# NOTE: Do NOT source ~/.bashrc or conda activate here.
# The Slurm script / interactive session must already have the conda environment
# activated before sourcing this config. Re-sourcing ~/.bashrc would re-initialize
# conda to `base` and undo the caller's conda activation.

# ============================================================================
# Environment
# ============================================================================
CONDA_ENV="bulk-HiC-processing"
R_CONDA_ENV="r-visualizations"

# ============================================================================
# Base directories
# ============================================================================
BASE_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
FIGURE_DIR="${BASE_DIR}/figure/cenpa-cuttag-enrichment"
CUTTAG_DIR="${BASE_DIR}/figure/cenpa-cuttag-centromere"

# ============================================================================
# Input files
# ============================================================================

# BAM files (k=1 unique mapping)
BAM_XG150="${CUTTAG_DIR}/mapped_reads/XG_150_aligned.sorted.bam"
BAM_XG151="${CUTTAG_DIR}/mapped_reads/XG_151_aligned.sorted.bam"
BAM_XG152="${CUTTAG_DIR}/mapped_reads/XG_152_aligned.sorted.bam"
BAM_XG153="${CUTTAG_DIR}/mapped_reads/XG_153_aligned.sorted.bam"

# BAM files (k=100 multimapping — for sensitivity analysis)
BAM_DIR_K100="${CUTTAG_DIR}/mapped_reads/k100-mapping"
BAM_XG150_K100="${BAM_DIR_K100}/XG_150_aligned.sorted.bam"
BAM_XG151_K100="${BAM_DIR_K100}/XG_151_aligned.sorted.bam"
BAM_XG152_K100="${BAM_DIR_K100}/XG_152_aligned.sorted.bam"
BAM_XG153_K100="${BAM_DIR_K100}/XG_153_aligned.sorted.bam"

# Sample metadata
declare -A SAMPLE_ANTIBODY
SAMPLE_ANTIBODY["XG_150"]="CENPA"
SAMPLE_ANTIBODY["XG_151"]="CENPA"
SAMPLE_ANTIBODY["XG_152"]="H3K27ac"
SAMPLE_ANTIBODY["XG_153"]="H3K27ac"

declare -A SAMPLE_REPLICATE
SAMPLE_REPLICATE["XG_150"]="rep1"
SAMPLE_REPLICATE["XG_151"]="rep2"
SAMPLE_REPLICATE["XG_152"]="control"
SAMPLE_REPLICATE["XG_153"]="control2"

# BAM list for iteration
BAM_LIST=("XG_150" "XG_151" "XG_152" "XG_153")
declare -A BAM_PATHS
BAM_PATHS["XG_150"]="$BAM_XG150"
BAM_PATHS["XG_151"]="$BAM_XG151"
BAM_PATHS["XG_152"]="$BAM_XG152"
BAM_PATHS["XG_153"]="$BAM_XG153"

# centroAnno output
CENTROANNO_DIR="${BASE_DIR}/output/outputs-from-centraAnno/hifiasm-0414/cautils-chrOnly"
REPEAT_REGIONS="${CENTROANNO_DIR}/repeat_regions.bed"
HORS_BED="${CENTROANNO_DIR}/HORs.bed"
CHR4_HORS="${CENTROANNO_DIR}/chr4_HORs.bed"
TOP10_REPEATS="${CENTROANNO_DIR}/top10_repeats.bed"
REPEAT_REGIONS_ALLSEQ="${BASE_DIR}/output/outputs-from-centraAnno/hifiasm-0414/cautils-all-seq/repeat_regions.bed"

# HiCAT chr4
HICAT_DIR="${BASE_DIR}/code/command-line-script/genome-annotation/HiCAT/300Gb"
HICAT_DECOMP="${HICAT_DIR}/final_decomposition.bed"
HICAT_TSV="${HICAT_DIR}/final_decomposition.tsv"

# Genome assembly
ASSEMBLY_DIR="${BASE_DIR}/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned-contamFiltered"
ASSEMBLY_FASTA="${ASSEMBLY_DIR}/assembly_final.sorted.headerRenamed.chrAssigned.contamFiltered.fasta"
ASSEMBLY_FAI="${ASSEMBLY_FASTA}.fai"

# JBAT gap annotations (for cross-reference)
JBAT_GAP_BED="${BASE_DIR}/output/outputs-from-haphic-alignment/references_hifiasm_male403_hifiHiCMode_041425_trail1_allChrom/04.build/out_JBAT_review_with_orig_name.bed"

# ============================================================================
# Output directories (created by scripts)
# ============================================================================
DATA_DIR="${FIGURE_DIR}/data"
FOREGROUND_DIR="${DATA_DIR}/foregrounds"
FRAGMENT_DIR="${DATA_DIR}/fragments"
EXCLUSION_DIR="${DATA_DIR}/exclusion"
BACKGROUND_DIR="${DATA_DIR}/backgrounds"
COUNTS_DIR="${DATA_DIR}/counts"
COVERAGE_DIR="${DATA_DIR}/coverage"
MAPPABILITY_DIR="${DATA_DIR}/mappability"
QC_DIR="${DATA_DIR}/qc"
DOMAIN_DIR="${DATA_DIR}/domains"
DOMAIN_FLANKS_DIR="${DATA_DIR}/domain_flanks"
DOMAIN_COUNTS_DIR="${DATA_DIR}/domain_counts"
DOMAIN_DISTANCE_DIR="${DATA_DIR}/domain_distance_profiles"
PLOTS_DIR="${FIGURE_DIR}/plots"
MAIN_PLOTS="${PLOTS_DIR}/main"
SUPP_PLOTS="${PLOTS_DIR}/supp"
RESULTS_DIR="${FIGURE_DIR}/results"

# ============================================================================
# Analysis parameters
# ============================================================================

# Chromosomes to include
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chrX chrY"

# Mapping parameters
MAPQ_THRESHOLD=10
SAMTOOLS_FILTER="-f 0x2 -F 0x900"   # proper pair, exclude secondary+supplementary

# Foreground Set B parameters (merge distance + flank padding)
# Primary
MERGE_DIST_1=50000
FLANK_PAD_1=10000
# Tighter
MERGE_DIST_2=25000
FLANK_PAD_2=5000
# Looser
MERGE_DIST_3=100000
FLANK_PAD_3=20000

# Domain definition parameters (V2)
DOMAIN_MERGE_DIST=250000       # merge CENP-A-positive intervals within 250 kb
DOMAIN_FLANK_DIST=1000000      # ±1 Mb flanks for local background
DOMAIN_MIN_INTERVALS=1         # minimum intervals per domain (1 = include singletons)

# BigWig files (for browser tracks in R)
BW_XG150="${CUTTAG_DIR}/bw_files/XG_150.all.bw"
BW_XG151="${CUTTAG_DIR}/bw_files/XG_151.all.bw"
BW_XG152="${CUTTAG_DIR}/bw_files/XG_152.all.bw"
BW_XG153="${CUTTAG_DIR}/bw_files/XG_153.all.bw"

# Background parameters
N_SHUFFLES=1000
SHUFFLE_SEED=20260716

# Local flank distances
FLANK_DISTANCES=(50000 100000 200000)

# Distance metaprofile bins (bp from interval boundary)
# Negative = upstream, 0 = within interval, Positive = downstream
DISTANCE_BINS=(
    "-100000:-50000"
    "-50000:-20000"
    "-20000:-10000"
    "-10000:0"
    "0"
    "0:10000"
    "10000:20000"
    "20000:50000"
    "50000:100000"
)

# ============================================================================
# Slurm settings (override in individual scripts if needed)
# ============================================================================
SLURM_PARTITION="condo"
SLURM_QOS="condo"
SLURM_ACCOUNT="csd788"
SLURM_CPUS=16
SLURM_MEM="32G"
SLURM_TIME="1-00:00:00"
SLURM_MAIL="you@example.com"
CLUSTER_LOG_DIR="/tscc/nfs/home/jhc103/cluster-logs"

# ============================================================================
# Helper functions
# ============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

check_file() {
    local f="$1"
    local label="${2:-file}"
    if [[ ! -f "$f" ]]; then
        log "ERROR: ${label} not found: $f"
        exit 1
    fi
    log "OK: ${label}: $f"
}

check_dir() {
    local d="$1"
    local label="${2:-directory}"
    if [[ ! -d "$d" ]]; then
        log "ERROR: ${label} not found: $d"
        exit 1
    fi
    log "OK: ${label}: $d"
}

# Create output directories
init_dirs() {
    mkdir -p "$FOREGROUND_DIR" "$FRAGMENT_DIR" "$EXCLUSION_DIR" \
             "$BACKGROUND_DIR/bg1_chrom_shuffle" \
             "$BACKGROUND_DIR/bg2_mappability_match" \
             "$BACKGROUND_DIR/bg3_local_flanks" \
             "$COUNTS_DIR" "$COVERAGE_DIR" "$MAPPABILITY_DIR" "$QC_DIR" \
             "$DOMAIN_DIR" "$DOMAIN_FLANKS_DIR" "$DOMAIN_COUNTS_DIR" \
             "$DOMAIN_DISTANCE_DIR" \
             "$MAIN_PLOTS" "$SUPP_PLOTS" "$RESULTS_DIR"
}

log "config.sh loaded successfully"
