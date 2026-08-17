#!/bin/bash
# config.sh — paths and parameters for the CENP-A domain-centered
# repeat-composition analysis (cenpa-repeat-chromosome).
# Source this from every script: source scripts/config.sh

set -euo

# ============================================================================
# Base directories
# ============================================================================
BASE_DIR="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
WORK_DIR="${BASE_DIR}/figure/cenpa-repeat-chromosome"
SRC_FIG_DIR="${BASE_DIR}/figure/cenpa-cuttag-enrichment"      # prior CENP-A analysis
CUTTAG_DIR="${BASE_DIR}/figure/cenpa-cuttag-centromere"       # BAMs / reads / BigWig

# Slurm / cluster
CLUSTER_LOG_DIR="/tscc/nfs/home/jhc103/cluster-logs"
SLURM_PARTITION="condo"
SLURM_QOS="condo"
SLURM_ACCOUNT="csd788"
SLURM_MAIL="you@example.com"

# Conda envs
CONDA_ENV="bulk-HiC-processing"      # samtools, bedtools, pysam
R_ENV="r-visualizations"

# ============================================================================
# Genome assembly (reference used for alignment; coordinates match BAM @SQ)
# ============================================================================
ASSEMBLY_DIR="${BASE_DIR}/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned-contamFiltered"
ASSEMBLY_FASTA="${ASSEMBLY_DIR}/assembly_final.sorted.headerRenamed.chrAssigned.contamFiltered.fasta"
ASSEMBLY_FAI="${ASSEMBLY_FASTA}.fai"

# ============================================================================
# CUT&Tag data
# ============================================================================
# k=1 fragment BEDs (chr start end, one per fragment)
FRAGMENTS_DIR="${SRC_FIG_DIR}/data/fragments"
FRAG_150="${FRAGMENTS_DIR}/XG_150_fragments.bed"   # CENP-A rep1
FRAG_151="${FRAGMENTS_DIR}/XG_151_fragments.bed"   # CENP-A rep2
FRAG_152="${FRAGMENTS_DIR}/XG_152_fragments.bed"   # H3K27ac ctrl
FRAG_153="${FRAGMENTS_DIR}/XG_153_fragments.bed"   # H3K27ac ctrl2

# k=100 BAMs (no NH tags; NH computed in step 1)
K100_DIR="${CUTTAG_DIR}/mapped_reads/k100-mapping"
BAM_150_K100="${K100_DIR}/XG_150_aligned.sorted.bam"
BAM_151_K100="${K100_DIR}/XG_151_aligned.sorted.bam"
BAM_152_K100="${K100_DIR}/XG_152_aligned.sorted.bam"

# trimmed reads (for k-mer leg)
READS_DIR="${CUTTAG_DIR}/input-reads"
READ_150_1="${READS_DIR}/XG_150_S1_L001_R1_001_val_1.fq.gz"
READ_150_2="${READS_DIR}/XG_150_S1_L001_R2_001_val_2.fq.gz"
READ_151_1="${READS_DIR}/XG_151_S2_L001_R1_001_val_1.fq.gz"
READ_151_2="${READS_DIR}/XG_151_S2_L001_R2_001_val_2.fq.gz"
READ_152_1="${READS_DIR}/XG_152_S3_L001_R1_001_val_1.fq.gz"
READ_152_2="${READS_DIR}/XG_152_S3_L001_R2_001_val_2.fq.gz"
READ_153_1="${READS_DIR}/XG_153_S1_L001_R1_001_val_1.fq.gz"
READ_153_2="${READS_DIR}/XG_153_S1_L001_R2_001_val_2.fq.gz"

declare -A SAMPLE_ANTIBODY=( [XG_150]=CENPA [XG_151]=CENPA [XG_152]=H3K27ac [XG_153]=H3K27ac )
declare -A SAMPLE_REP=( [XG_150]=rep1 [XG_151]=rep2 [XG_152]=ctrl [XG_153]=ctrl2 )

# ============================================================================
# Repeat arrays (TRF period bins, merged)
# ============================================================================
PERIOD_DIR="${SRC_FIG_DIR}/period-enrichment"
ARRAYS_DIR="${PERIOD_DIR}/data/merged/arrays"
BIN4_195="${ARRAYS_DIR}/bin4_arrays.bed"   # 193-195 bp (L1 5' tandem)  chr start end id cnt
BIN6_349="${ARRAYS_DIR}/bin6_arrays.bed"   # 348-349 bp (satellite)
BIN8_389="${ARRAYS_DIR}/bin8_arrays.bed"   # 386-390 bp (195 dimer)

# All TRF loci (period bins) for repeat-density track
TRF_BED="${PERIOD_DIR}/data/trf_chr_only_period_bins.bed"

# L1 family overlap for 195 arrays (RepeatMasker; seq coords -> chr map)
RM_OVERLAP_DIR="${SRC_FIG_DIR}/195bp-repeatmasker-overlap"
RM_195_FAMILY="${RM_OVERLAP_DIR}/bin4_vs_rm_family_wo.tsv"
SEQ2CHR="${RM_OVERLAP_DIR}/seq2chr_map.tsv"

# Consensus sequences for k-mer probe definitions
L1_189_CONS="${RM_OVERLAP_DIR}/rnd-1_family-189.consensus.fa"
L1_18_CONS="${RM_OVERLAP_DIR}/rnd-1_family-18.consensus.fa"
TANDEM_349_CONS="${BASE_DIR}/code/command-line-script/tandem-repeat-investigation/349peak_repeat_10Longest_consensusSeq.fasta"
TANDEM_195_CONS="${BASE_DIR}/code/command-line-script/tandem-repeat-investigation/195peak_repeat_10Longest_consensusSeq.fasta"

# ============================================================================
# Chromosomes
# ============================================================================
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chrX chrY"
CHROM_SIZES="${SRC_FIG_DIR}/data/chrom_sizes.txt"

# ============================================================================
# Analysis parameters
# ============================================================================
WINDOW=100000                 # 100 kb non-overlapping windows
BAND_NEAR=250000             # core -> 250 kb
BAND_FAR=1000000             # 250 kb -> 1 Mb
K31=31                        # k-mer length for probes
N_SHUFFLE=200                 # matched-null placements per array
SHUFFLE_SEED=20260807

# ============================================================================
# Output directories (created by init_dirs)
# ============================================================================
DATA_DIR="${WORK_DIR}/data"
NH_DIR="${DATA_DIR}/nh"
COVW_DIR="${DATA_DIR}/coverage_weighted"
DOMAIN_DIR="${DATA_DIR}/domains"
MAPPABILITY_DIR="${DATA_DIR}/mappability"
GC_DIR="${DATA_DIR}/gc"
REPDENS_DIR="${DATA_DIR}/repdensity"
PROBE_DIR="${DATA_DIR}/probes"
RESULTS_DIR="${WORK_DIR}/results"
MAIN_PLOTS="${WORK_DIR}/plots/main"
SUPP_PLOTS="${WORK_DIR}/plots/supp"
SCRIPTS_DIR="${WORK_DIR}/scripts"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

check_file() { local f="$1"; local label="${2:-file}"; if [[ ! -f "$f" ]]; then log "ERROR: ${label} not found: $f"; exit 1; fi; log "OK: ${label}: $f"; }

init_dirs() {
    mkdir -p "$DATA_DIR" "$NH_DIR" "$COVW_DIR" "$DOMAIN_DIR" "$MAPPABILITY_DIR" \
             "$GC_DIR" "$REPDENS_DIR" "$PROBE_DIR" "$RESULTS_DIR" \
             "$MAIN_PLOTS" "$SUPP_PLOTS" "$SCRIPTS_DIR"
}

log "config.sh loaded (work dir: $WORK_DIR)"
