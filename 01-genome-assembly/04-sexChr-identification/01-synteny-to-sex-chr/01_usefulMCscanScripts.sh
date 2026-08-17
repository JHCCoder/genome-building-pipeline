#!/bin/bash
# Synteny to sex chromosomes using JCVI / MCScanX.
#
# This is a runbook (run the commands in order), not a single monolithic job.
# Idea: transfer annotation to the de-novo assembly, then find synteny between
# our assembly and a reference with sex-chromosome information (here mouse
# GRCm39). Scaffolds with strong/complete synteny to an annotated sex chromosome
# are sex-linked.
#
# Required helper scripts (expected in this directory, not yet in the repo):
#   count_transcript.sh  make_seqid.sh  make_layout.sh  highlight_chrom.sh
#   highlight_chrom1.sh  find_black.py
# Required data files (expected in the working directory):
#   GCF_000001635.27_GRCm39_genomic.gff  chrom_dict_mouse_GRC29.txt  blocks.layout

# Load shared configuration (config.sh at the repo root)
_repo_root="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [[ ! -f "$_repo_root/config.sh" && "$_repo_root" != "/" ]]; do _repo_root="$(dirname "$_repo_root")"; done
source "$_repo_root/config.sh"
unset _repo_root

source ~/.bashrc   # initialize conda (adjust for your setup)
conda activate "$ENV_JCVI"

## Process mouse data
gffread GCF_000001635.27_GRCm39_genomic.gff -g "$MOUSE_GENOME_FA" -x GRCm39.cds.fna
python -m jcvi.formats.fasta format GRCm39.cds.fna mouse.cds
python -m jcvi.formats.gff bed --type=mRNA --key=ID input/GCF_000001635.27_GRCm39_genomic.gff -o mouse_orig.bed
# Rename chromosomes using a dict built from the NCBI page for the genome
awk -F'\t' -v OFS='\t' 'NR==FNR {dict[$2]=$1; next} $1 in dict {$1=dict[$1]} 1' chrom_dict_mouse_GRC29.txt mouse_orig.bed > mouse.bed

## Process degu data
gffread "$DEGU_LIFTOFF_GFF" -g "$DEGU_SCAFFOLDS_FA" -x hifiasm_041425_haphic.cds.fna
python -m jcvi.formats.fasta format hifiasm_041425_haphic.cds.fna degus.cds
python -m jcvi.formats.gff bed --type=mRNA --key=ID "$DEGU_LIFTOFF_GFF" -o degus.bed

## Run synteny, dotplot, and per-chromosome gene depth
python -m jcvi.compara.catalog ortholog mouse degus --cscore=.99 --no_strip_names
python -m jcvi.graphics.dotplot mouse.degus.anchors
python -m jcvi.compara.synteny depth --histogram mouse.degus.anchors

## Check whether an anchor exists for a particular chromosome
./count_transcript.sh ChrY mouse.bed mouse.degus.anchors

#### MACROSYNTENY (needs seqids + layout files) ####
python -m jcvi.compara.synteny screen --minspan=30 --simple mouse.degus.anchors mouse.degus.anchors.simple
./make_seqid.sh mouse.bed degus.bed mouse.degus.anchors.simple seqids
./make_layout.sh mouse.bed degus.bed mouse.degus.anchors.simple layout
python -m jcvi.graphics.karyotype seqids layout
python -m jcvi.graphics.karyotype seqids layout --figsize 14x10 --dpi 300 --o large_karyotype --format png

## Highlight certain alignments
./highlight_chrom.sh mouse.degus.anchors.simple mouse.bed degus.bed mouse.degus.anchors.chrX.simple "ChrX:r"
./make_layout.sh mouse.bed degus.bed mouse.degus.anchors.chrX.simple layout
python -m jcvi.graphics.karyotype seqids layout --figsize 14x10 --dpi 300 --o karyotype_chrX.pdf --chrstyle=roundrect

#### MICROSYNTENY ####
python -m jcvi.compara.synteny mcscan mouse.bed mouse.degus.lifted.anchors --iter=1 -o mouse.degus.i1.blocks
python find_black.py input/GCF_000001635.27_GRCm39_genomic.gff mouse.degus.i1.blocks App app_blocks.txt
python -m jcvi.formats.bed merge mouse.bed degus.bed -o mouse_degus.bed
# cp ../hifi022425scaf-GRCm38-mouse-degu/blocks.layout .   # context-specific layout reuse
python -m jcvi.graphics.synteny app_blocks.txt mouse_degus.bed blocks.layout --glyphstyle=arrow --format png
