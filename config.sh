#!/usr/bin/env bash
# =============================================================================
# Shared configuration for the genome-building pipeline.
#
# This file is sourced by every pipeline script. Edit the paths below to match
# your own environment and cluster; the pipeline scripts themselves should not
# need to be changed.
# =============================================================================

# --- Scheduler (Slurm on TSCC/UCSD — adjust for your cluster) ----------------
export EMAIL="you@example.com"              # your email for job notifications

# --- Top-level directories ----------------------------------------------------
export PROJ_ROOT="/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
export DATA_DIR="$PROJ_ROOT/data"
export OUTPUT_DIR="$PROJ_ROOT/output"
export CODE_DIR="$PROJ_ROOT/code"
export TOOLSHED_DIR="/tscc/projects/ps-renlab2/jhc103/toolshed"   # locally installed tools
export SCRATCH_DIR="/tscc/lustre/ddn/scratch/jhc103"              # scratch / temp space
export CLUSTER_LOG_DIR="/tscc/nfs/home/jhc103/cluster-logs"       # Slurm %x.%j.%N logs

# --- Conda environments -------------------------------------------------------
export ENV_GENOME_ANNOTATION="genome-annotation"
export ENV_BULK_HIC="bulk-HiC-processing"
export ENV_HAPHIC="haphic"
export ENV_DISCOVERY="toolshed-DiscoverY"
export ENV_JCVI="toolshed-jcvi"

# --- Tool binaries (invoked by absolute path) --------------------------------
export HIFIASM_BIN="$TOOLSHED_DIR/hifiasm/hifiasm"
export HAPHIC_DIR="$TOOLSHED_DIR/HapHiC"           # contains 'haphic' and 'utils/filter_bam'

# =============================================================================
# 01-genome-assembly — stage-specific inputs
# =============================================================================

# HiFi reads (hifiasm)
export HIFI_READS="$SCRATCH_DIR/hifi/male/*/*fastq"

# Hi-C reads — raw trim_galore ("_val_") output, used by hifiasm --h1/--h2
export HIC_RAW_R1="$DATA_DIR/sequencing-reads-HiC/WB_438_2_S1_L008_R1_001_val_1.fq.gz"
export HIC_RAW_R2="$DATA_DIR/sequencing-reads-HiC/WB_438_2_S1_L008_R2_001_val_2.fq.gz"

# Hi-C reads — combined/trimmed, used by HapHiC
export HIC_HAPHIC_R1="$CODE_DIR/command-line-script/haphic/input-hic-read/WB_438_R1_trimmed_combined.fastq.gz"
export HIC_HAPHIC_R2="$CODE_DIR/command-line-script/haphic/input-hic-read/WB_438_R2_trimmed_combined.fastq.gz"

# hifiasm output
export HIFIASM_OUT_DIR="$OUTPUT_DIR/outputs-from-hifiasm"
export HIFIASM_ASM_NAME="deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_052325.asm"

# Reference assembly to scaffold with HapHiC
export HAPHIC_REF="$DATA_DIR/OctDegus1_genome/OctDeg1/fasta/genome.fa"

# Bulk Hi-C preprocessing (bwa-mem2 + pairtools + cooler) inputs/outputs
export BHIC_GENOME="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta"
export BHIC_GENOME_FAI="${BHIC_GENOME}.fai"
export BHIC_RAW_INPUT_DIR="/tscc/projects/ps-renlab2/nzemke/Element_share/20250811"   # raw Hi-C fastq source
export BHIC_DATA_DIR="$DATA_DIR/sequencing-reads-HiC"
export BHIC_MAP_DIR="$OUTPUT_DIR/mapped_alignments"
export BHIC_MATRIX_DIR="$OUTPUT_DIR/matrix"
export BHIC_FILE_LIST="$DATA_DIR/file_list_111824_deepSeq_403Male.txt"

# Mitochondria filtering
export MITO_FASTA="$DATA_DIR/denovo_OctDegus_genome/degu_mitochondria/octDeg1_mitochondria_chromosome.fasta"
export MITO_BLAST_ASSEMBLY="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly/deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_041325.asm.hic.p_ctg.fa"

# Sex-chromosome identification
export SRY_GENE="$CODE_DIR/command-line-script/assigning-chromosomes/sry-gene/all_sry_gene.fasta"
export SRY_BLAST_ASSEMBLY="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked/assembly_final.sorted.headerRenamed.fasta.masked"

# Synteny to sex chromosomes (JCVI / MCScanX)
export MOUSE_GENOME_FA="$DATA_DIR/GRCm39_genome/GCF_000001635.27_GRCm39_genomic.fna"
export DEGU_LIFTOFF_GFF="$OUTPUT_DIR/outputs-from-liftoff/hifiasm-041425-scaffolded/hifiasm-041425-scaffolded.gff"
export DEGU_SCAFFOLDS_FA="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded/scaffolds.fa"

# =============================================================================
# 02-initial-genome-assembly-evaluation — Merqury + BUSCO inputs
# =============================================================================
export ENV_MERQURY="toolshed-merqury"
export ENV_BUSCO="toolshed-busco"

# Final assembly to evaluate (contam-filtered + haplotig-purged)
export FINAL_ASSEMBLY="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned-contamFiltered-purgeDupSeq-haplotigOnly/assembly_final.sorted.headerRenamed.chrAssigned.contamFiltered.purgeDupSeqHaplotigOnly.fasta"

# Merqury (k-mer-based completeness/consensus)
export MERQURY_OUT_DIR="$OUTPUT_DIR/outputs-from-merqury"
export MERYL_DB="collaborator_degu_WGS.meryl"     # Illumina WGS k-mer db (see 01_meryl_db.sh)
export WGS_READS="/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/data/sequencing-illumina/from-WGS-collaborator/*/*.fastq.gz"

# BUSCO (gene-content completeness)
export BUSCO_OUT_DIR="$OUTPUT_DIR/outputs-from-busco-ortholog-alignment"
export BUSCO_LINEAGE="glires_odb10"
