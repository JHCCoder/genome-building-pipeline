#!/bin/bash

source ~/.bashrc
conda activate toolshed-jcvi 


## Some gffs for degus:
# 1. braker gff: ~/ps-renlab2-link/degu-genome-assembly-proj/code/command-line-script/genome-annotation/annotate-braker3-results/annotate-gff/braker.uniprotBlast.interproscan.gff/braker.gff3


## Process mouse data
gffread GCF_000001635.27_GRCm39_genomic.gff -g /tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/GRCm39_genome/GCF_000001635.27_GRCm39_genomic.fna -x GRCm39.cds.fna


python -m jcvi.formats.fasta format GRCm39.cds.fna mouse.cds

python -m jcvi.formats.gff bed --type=mRNA --key=ID input/GCF_000001635.27_GRCm39_genomic.gff -o mouse_orig.bed

#awk 'NR==FNR {dict[$2]=$1; next} $1 in dict {$1=dict[$1]} 1' chrom_dict_mouse_GRC29.txt mouse_orig.bed > mouse.bed # create the dict using the ncbi page for the genome 
awk -F'\t' -v OFS='\t' 'NR==FNR {dict[$2]=$1; next} $1 in dict {$1=dict[$1]} 1' chrom_dict_mouse_GRC29.txt mouse_orig.bed > mouse.bed

## Process degu data 
gffread /tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/output/outputs-from-liftoff/hifiasm-041425-scaffolded/hifiasm-041425-scaffolded.gff -g ~/ps-renlab2-link/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded1/scaffolds.fa -x hifiasm_041425_haphic.cds.fna

python -m jcvi.formats.fasta format hifiasm_041425_haphic.cds.fna degus.cds

python -m jcvi.formats.gff bed --type=mRNA --key=ID /tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/output/outputs-from-liftoff/hifiasm-041425-scaffolded/hifiasm-041425-scaffolded.gff -o degus.bed

## Run the synteny
python -m jcvi.compara.catalog ortholog mouse degus --cscore=.99 --no_strip_names

## Run the dotplot 
python -m jcvi.graphics.dotplot mouse.degus.anchors

## Run the depth or the # of genes match
python -m jcvi.compara.synteny depth --histogram mouse.degus.anchors

##  Checking whether the anchor is available for a particular chromosome
./count_transcript.sh ChrY mouse.bed mouse.degus.anchors 

#### Examine MACROSYNTENY ######
#### Require 2 additional files: seqids => sets of chromosomes to include
####				 layout => how to draw the chromosomes
## Make more simple anchors 
python -m jcvi.compara.synteny screen --minspan=30 --simple mouse.degus.anchors mouse.degus.anchors.simple

## Make a karyotype seqids file first
./make_seqid.sh mouse.bed degus.bed mouse.degus.anchors.simple seqids

## Make a layout file 
./make_layout.sh mouse.bed degus.bed mouse.degus.anchors.simple layout

## Write synteny plot
python -m jcvi.graphics.karyotype seqids layout
## You can also add parameter to command:
python -m jcvi.graphics.karyotype seqids layout --figsize 14x10 --dpi 300 --o large_karyotype --format png

## ==> Highlight certain alignments 
./highlight_chrom.sh mouse.degus.anchors.simple mouse.bed degus.bed mouse.degus.anchors.chrX.simple "ChrX:r"
./make_layout.sh mouse.bed degus.bed mouse.degus.anchors.chrX.simple layout
python -m jcvi.graphics.karyotype seqids layout --figsize 14x10 --dpi 300 --o karyotype_chrX.pdf --chrstyle=roundrect
## Or multiple alignments
./highlight_chrom1.sh mouse.degus.anchors.simple mouse.bed degus.bed mouse.degus.anchors.chr1.chrX.simple1 "Chr1:g*" "ChrX:r*"
# Then change your layout file and you rerun the karyotype command


#### Examine Microsynteny ####

## Isolate the regions 
python -m jcvi.compara.synteny mcscan mouse.bed mouse.degus.lifted.anchors --iter=1 -o mouse.degus.i1.blocks

#App => rna-XR_380384.5, rna-XM_006522873.4 *, rna-NM_007471.3, rna-NM_001198825.1, rna-NM_001198824.1, rna-NM_001198823.1, rna-XM_006522874.3, rna-XM_036159749.1
python find_black.py input/GCF_000001635.27_GRCm39_genomic.gff mouse.degus.i1.blocks App app_blocks.txt

python -m jcvi.formats.bed merge mouse.bed degus.bed -o mouse_degus.bed

cp ../hifi022425scaf-GRCm38-mouse-degu/blocks.layout .

python -m jcvi.graphics.synteny app_blocks.txt mouse_degus.bed blocks.layout --glyphstyle=arrow --format png

