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
export HIFI_READS="$SCRATCH_DIR/hifi/male/*/*fastq.gz"

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

# Read-coverage mapping (bwa-mem2 short-read / minimap2 long-read -> bedtools coverage)
export ENV_BEDTOOLS="toolshed-bedtools"
export COVERAGE_ASSEMBLY="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta"
export COVERAGE_SHORT_R1="$DATA_DIR/sequencing-illumina/from-WGS-collaborator/SRR19145623_24threads/SRR19145623_1_val_1.fq.gz $DATA_DIR/sequencing-illumina/from-WGS-collaborator/SRR19145629_24threads/SRR19145629_1_val_1.fq.gz"
export COVERAGE_SHORT_R2="$DATA_DIR/sequencing-illumina/from-WGS-collaborator/SRR19145623_24threads/SRR19145623_2_val_2.fq.gz $DATA_DIR/sequencing-illumina/from-WGS-collaborator/SRR19145629_24threads/SRR19145629_2_val_2.fq.gz"
export COVERAGE_HIFI_READS="$SCRATCH_DIR/hifi/male/*/*fastq.gz"
export COVERAGE_OUT_DIR="$OUTPUT_DIR/outputs-from-read-coverage"

# GenomeScope (KMC k-mer histogram -> genome size / heterozygosity / repeat estimate)
export KMC_BIN="$TOOLSHED_DIR/kmc/bin/kmc"
export KMC_TOOLS_BIN="$TOOLSHED_DIR/kmc/bin/kmc_tools"
export GENOMESCOPE_READS="$DATA_DIR/sequencing-illumina/male/processed-reads/*fastq.gz"
export GENOMESCOPE_OUT_DIR="$OUTPUT_DIR/outputs-from-genomeScope-kmc"
export GENOMESCOPE_TEMP_DIR="$SCRATCH_DIR/temp-dir-fast"

# =============================================================================
# 03-genome-annotation — repeat masking + annotation inputs
# =============================================================================
export ENV_REPEATMODELER="toolshed-repeatmodeler"
export ENV_LIFTOFF="toolshed-liftoff"

# Final masked / mito-assigned assemblies (annotation target)
export MASKED_ASSEMBLY="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.masked.fasta"
export MITO_ASSEMBLY="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta"

# Repeat masking (RepeatModeler / RepeatMasker)
export REPEAT_MODELER_ASM="$OUTPUT_DIR/outputs-from-haphic-alignment/references_hifiasm_male403_hifiHiCMode_022425/04.build/scaffolds.fa"
export REPEAT_MASKER_ASM="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-hap-contamFiltered-mitoFiltered-scaffolded/hap2_scaffolds.fa"
export REPEAT_OUT_DIR="$OUTPUT_DIR/outputs-from-repeatmasker"

# Liftoff (transfer-based annotation)
export LIFTOFF_QUERY="$OUTPUT_DIR/outputs-from-repeatmasker/verkko_hap1_scaffolded/scaffolds.fa"
export LIFTOFF_REF="$MITO_ASSEMBLY"
export LIFTOFF_REF_GFF="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm_041425_denovoEnhanced_peaks2utr_sorted.gff3"

# De-novo (Braker3 + RNA-seq) inputs
export MRNA_DIR="$DATA_DIR/sequencing-mRNAseq/SRA-ncbi-deposits"
export HISAT2_OUT_DIR="$OUTPUT_DIR/outputs-from-hisat2Aligned-samtoolSorted-mRNA-transcript"
export HISAT2_INDEX_NAME="hifiasm_041425_haphic_masked_curated"
export BRAKER3_DIR="$CODE_DIR/command-line-script/genome-annotation/braker3"
export BRAKER3_SIF="$BRAKER3_DIR/braker3.sif"
export BRAKER3_PROTEINS="$BRAKER3_DIR/Vertebrata_plus_relatives.fa"

# Functional annotation tools
export BLASTP_BIN="$TOOLSHED_DIR/blast+/ncbi-blast-2.16.0+/bin/blastp"
export BLAST_DB_DIR="/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/data/protein-databases/uniprot-vertebrate"
export INTERPROSCAN_BIN="$TOOLSHED_DIR/interproscan/interproscan-5.74-105.0/interproscan.sh"

# =============================================================================
# 06-annotation-enhancement-peaks2utr-agat — 3' UTR annotation + finalization
# =============================================================================
export ENV_PEAKS2UTR="toolshed-peaks2utr"
export ENV_AGAT="toolshed-agat"
export ENV_DEEPTOOLS="toolshed-deeptools"        # bamCoverage (expression BigWig)

# scMultiome GEX data (4 samples merged: 181_PFC, 181_dHIP, 6997_dHIP, 7000_PFC)
export PEAKS2UTR_BAM="$CODE_DIR/command-line-script/genome-annotation/peaks2utr/merged.bam"
export PEAKS2UTR_BW="$CODE_DIR/command-line-script/genome-annotation/peaks2utr/merged.bw"
export PEAKS2UTR_WORK_DIR="$CODE_DIR/command-line-script/genome-annotation/peaks2utr"

# Annotation files for the UTR + AGAT finalization
export MERGE_OUT_DIR="$CODE_DIR/command-line-script/annotation-merging/output"
export MERGED_GFF="$MERGE_OUT_DIR/hifiasm-041425-denovoEnhanced_merged.gff3"           # 05 merge output
export PEAKS2UTR_PREV_INPUT="$PEAKS2UTR_WORK_DIR/hifiasm_041425_denovoEnhanced_sorted.gff"   # previous annotation (full peaks2utr input)
export PEAKS2UTR_PREV_OUTPUT="$PEAKS2UTR_WORK_DIR/hifiasm_041425_denovoEnhanced_peaks2utr_sorted_perfect.gff"  # previous full peaks2utr output (AGAT-processed)
export CHANGED_SUBSET_GFF="$MERGE_OUT_DIR/peaks2utr_changed_subset.gff3"              # changed-gene subset (02)
export PEAKS2UTR_SUBSET_OUT="$MERGE_OUT_DIR/peaks2utr_subset_out.gff3"                # subset peaks2utr output (03)
export FINAL_UTR_GFF="$MERGE_OUT_DIR/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.gff3"   # assembled w/ UTRs (04)
export FINAL_AGAT_GFF="$MERGE_OUT_DIR/hifiasm-041425-denovoEnhanced_peaks2utr_sorted.agat.gff3"  # final AGAT-normalized (05)

# =============================================================================
# 05-genome-assembly-indepth-assessment-contam-filtering — purge_dup inputs
# =============================================================================
export ENV_ASSESSMENT="genome-assembly-assessment"
export PURGE_DUPS_BIN="$TOOLSHED_DIR/purge_dups/bin"

# Scaffolded assembly to purge (chromosome-assigned, unmasked) + its haplotig assembly
export CHR_ASSIGNED_ASSEMBLY="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.fasta"
export PURGE_HAP_ASM="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly/deNovo_hifiHiCMode_hifiData_aggressivePurge3_kmer21_052325.asm.hic.a_ctg.fa"
export PURGE_OUT_DIR="$OUTPUT_DIR/outputs-from-purge-duplicate"
export PURGE_WORK_DIR="$PURGE_OUT_DIR/hifiasm-041425-scaffolded-assembly"
export PURGE_PB_LIST="$SCRATCH_DIR/hifi/male/long-reads/m54284U_220310_223856.hifi_reads.fastq.gz $SCRATCH_DIR/hifi/male/long-reads/m54284U_230214_161608.hifi_reads.fastq.gz $SCRATCH_DIR/hifi/male/long-reads/m54284U_230217_165146.hifi_reads.fastq.gz $SCRATCH_DIR/hifi/male/long-reads/m84216_250214_011813_s2.fastq.gz"

# Final assembly + BISER segdup output (for the segdup enrichment analysis)
export DUP_BED="$PURGE_WORK_DIR/dups.bed"
export BISER_BEDPE="$CODE_DIR/command-line-script/genome-annotation/biser/hifiasm-041425/segdup_output_mod.bedpe"
export ASSEMBLY_MITO_FA="$PROJ_ROOT/degu-genome-browser-pythonVersion/assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta"

# =============================================================================
# 05-… — paralog read-depth screen inputs
# =============================================================================
export ENV_TRANSCRIPTOME_MAPPING="transcriptome-mapping"
export READ_DEPTH_ASSEMBLY="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned-contamFiltered/assembly_final.sorted.headerRenamed.chrAssigned.contamFiltered.fasta"
export READ_DEPTH_HIFI_DIR="$SCRATCH_DIR/hifi/male/long-reads"
export PARALOG_GFF="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm_041425_denovoEnhanced_peaks2utr_sorted_perfect.gff"
export PARALOG_REPEAT_BED="$PROJ_ROOT/figure/circos-plot/feature-overview/assembly_final.sorted.headerRenamed.fasta.out.gff"
export PARALOG_GENOME_LENGTHS="$CODE_DIR/command-line-script/contig-coverage/hifiasm_041425_scaffolded_juiceBox_sorted_hardMasked_chrAssigned_mitoAdded_genome_length.txt"
export STAR_INDEX="$PROJ_ROOT/figure/paralog-alignment-visualization/star_index"

# Genome polishing (meryl + winnowmap + yak + NextPolish2)
export ENV_GENOME_POLISHING="genome-polishing"

# Assembly to polish (mito + contam filtered, pre-scaffolding) and its reads
export POLISH_ASM="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-contamFiltered/genome_chrom_contamRemoved.fasta"
export POLISH_HIFI_READS="$SCRATCH_DIR/hifi/male/long-reads/*.fastq.gz"
export POLISH_ILLUMINA_R1="$DATA_DIR/sequencing-illumina/male/processed-reads/nR184-L3-G3-P20-GATTCCTT-TGCCTATG-READ1-Sequences-trimmed.fastq.gz"
export POLISH_ILLUMINA_R2="$DATA_DIR/sequencing-illumina/male/processed-reads/nR184-L3-G3-P20-GATTCCTT-TGCCTATG-READ2-Sequences-trimmed.fastq.gz"
export POLISH_WORK_DIR="$OUTPUT_DIR/outputs-from-nextpolish2"

# =============================================================================
# 07-segmental-duplication-analysis — BISER inputs
# =============================================================================
# BISER (user-level pip install); ENV_GENOME_ANNOTATION supplies samtools
export BISER_BIN="$HOME/.local/bin/biser"
export BISER_OUT_DIR="$OUTPUT_DIR/outputs-from-biser"
export BISER_THREADS=16
export BISER_GC_HEAP="5G"
export BISER_MAX_EDIT_ERROR=10
export BISER_MAX_ERROR=15
export BISER_GENOME_DIR="$DATA_DIR/genome-related-species"      # multi-species assemblies
export BISER_GENOME_LIST="$BISER_GENOME_DIR/genome_list7.txt"   # <assembly.fna>\t<species>

# =============================================================================
# 06-tandem-repeat-analysis — processing tools
# =============================================================================
export TRF_BIN="$TOOLSHED_DIR/TRF-4.09.1/build/src/trf"
export CENTROANNO_DIR="$TOOLSHED_DIR/centroAnno"
export HICAT_BIN="/tscc/projects/ps-renlab2/jhc103/miniconda3-storage/envs/toolshed-HiCAT/bin/hicat"

# =============================================================================
# 04-… — synteny (ntSynt) + mashtree inputs
# =============================================================================
export ENV_PHYLO="evolutionary-tree"            # perl-BioPerl (mashtree) + FastME

# mash / mashtree (Jaccard-distance tree ordering the synteny plot)
export MASH_BIN="$TOOLSHED_DIR/mash-Linux64-v2.3/mash"
export MASHTREE_DIR="$TOOLSHED_DIR/mash-Linux64-v2.3/mashtree"   # bin/{mashtree,mashtree_jackknife.pl,mashtree_bootstrap.pl}

# ntSynt (synteny blocks) + ntSynt-viz (ribbon plot); both run in ENV_ASSESSMENT
export NTSYNT_BIN="$TOOLSHED_DIR/ntSynt/bin/ntSynt"
export NTSYNT_VIZ_BIN="$TOOLSHED_DIR/ntsynt-viz/ntSynt-viz-1.0.0/bin/ntsynt_viz.py"

export SYNTENY_OUT_DIR="$OUTPUT_DIR/outputs-from-synteny"
export SYNTENY_FASTA_LIST="$CODE_DIR/github-code-to-share/04-genome-annotation-initial-evaluation/06-mashTree-syntenic-evaluation/fasta_list.txt"
export SYNTENY_NAME_CONVERSIONS="$CODE_DIR/github-code-to-share/04-genome-annotation-initial-evaluation/06-mashTree-syntenic-evaluation/name-conversions.tsv"
export SYNTENY_TARGET_GENOME="Degus"

# =============================================================================
# 04-genome-annotation-initial-evaluation — evolutionary-placement inputs
# =============================================================================
# BUSCO (genome mode, --metaeuk) runs reuse ENV_BUSCO (toolshed-busco) and
# ENV_BUSCO_LINEAGE defined below. The tree itself is built with MAFFT + trimAl
# + RAxML in the evolutionary-tree env.
export ENV_EVOLUTIONARY_TREE="evolutionary-tree"   # /tscc/.../miniconda3-storage/envs/evolutionary-tree (MAFFT 7.526, RAxML 8.2.12)
export EVO_OUT_DIR="$OUTPUT_DIR/outputs-from-evolutionary-placement"
export EVO_BUSCO_LINEAGE="euarchontoglires_odb12"  # BUSCO lineage for the species tree (auto-downloaded on first run)
export EVO_SPECIES_LIST="$CODE_DIR/github-code-to-share/04-genome-annotation-initial-evaluation/05-evolutionary-placement/species_list.tsv"

# =============================================================================
# 02-genome-assembly-initial-evaluation — Hi-C contact-map visualization
# =============================================================================
export ENV_HICEXPLORER="toolshed-HiCExplorer"   # hicBuildMatrix / hicPlotMatrix / hicPCA / hicMerge*
export ENV_SAMTOOLS="toolshed-samtools"         # samtools + sambamba

# Genome + restriction sites used for Hi-C matrix building (HiCExplorer)
export HIC_GENOME="$DATA_DIR/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned/assembly_final.sorted.headerRenamed.chrAssigned.fasta"
export HIC_RESTRICTION_BED="$PROJ_ROOT/figure/hic-plot/041425_assembly_GATC.bed"   # hicFindRestSite --fasta $HIC_GENOME --searchPattern GATC -o ...
export HIC_OUT_DIR="$PROJ_ROOT/figure/hic-plot"                                    # working/output dir (matrices, plots)
export HIC_HIFI_READ_DIR="$CODE_DIR/command-line-script/haphic/input-hic-read"    # HiFi Hi-C fastq
export HIC_SHORTREAD_DIR="$DATA_DIR/sequencing-illumina/from-WGS-collaborator"    # Illumina Hi-C fastq

# =============================================================================
# 05-… — structural-error detection (Inspector) inputs
# =============================================================================
export ENV_INSPECTOR="toolshed-inspector"          # inspector.py + minimap2 + samtools
export INSPECTOR_BIN="$TOOLSHED_DIR/Inspector/inspector.py"
export INSPECTOR_ASM="$MITO_ASSEMBLY"              # scaffolded, chr-assigned, mito-added assembly
export INSPECTOR_HIFI_READS="$SCRATCH_DIR/hifi/male/*/*fastq.gz"   # unquoted in the script to expand the glob
export INSPECTOR_REF="$DATA_DIR/OctDegus1_genome/OctDeg1/fasta/genome.fa"  # optional reference (reference-based mode)
export INSPECTOR_OUT_DIR="$OUTPUT_DIR/outputs-from-inspector"
