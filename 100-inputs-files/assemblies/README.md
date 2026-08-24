# assemblies/ — de novo assemblies (symlinks only)

The entries here are **symbolic links** to assembly FASTA files on the TSCC
cluster. Each is roughly 3 GB, so none of it is stored in this repository.

These are outputs of this study's pipeline, built from the raw reads rather than
downloaded. To regenerate them, obtain the reads and run stages 01–05.

| Data | Archive | Accession |
|---|---|---|
| HiFi long-read + Illumina short-read (assembly inputs) | NCBI BioProject | **PRJNA1393203** |
| Hi-C (scaffolding input) | NCBI GEO | **GSE314503** |

- NCBI BioProject: https://www.ncbi.nlm.nih.gov/bioproject/
- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/

While these accessions remain under embargo, reviewer links and GEO access
tokens are issued with the manuscript. Those credentials are **not stored in
this repository**.

## Contents

| Link name | What | `config.sh` var |
|---|---|---|
| `assembly_final.chrAssigned.mito.fasta` | chromosome-assigned + mitochondrion | `MITO_ASSEMBLY` |
| `assembly_final.chrAssigned.fasta` | chromosome-assigned | `CHR_ASSIGNED_ASSEMBLY` |
| `assembly_final.chrAssigned.masked.fasta` | repeat-masked | `MASKED_ASSEMBLY` |
| `assembly_final.chrAssigned.contamFiltered.fasta` | contamination-filtered | `READ_DEPTH_ASSEMBLY` |
| `assembly_final.chrAssigned.contamFiltered.purgeDupSeqHaplotigOnly.fasta` | + haplotig-purged | `FINAL_ASSEMBLY` |
| `hap2_scaffolds.fa` | hap2 scaffolds | — |
| `verkko_hap1_scaffolds.fa` | verkko hap1 scaffolds | `VERKKO_HAP1_SCAFFOLDS_FA` |
| `haphic_assembly_final.fasta` | HapHiC scaffolding build | `HAPHIC_ASSEMBLY_FA` |

An interactive genome browser for the OctDeg2.0 assembly is available at:

https://renlab.sdsc.edu/share/jhc103/degu-genome-browser/
