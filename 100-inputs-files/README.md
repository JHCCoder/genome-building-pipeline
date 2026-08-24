# Inputs (stand-in)

This directory is a **stand-in** for the input data used by the pipeline. The
entries are **symbolic links** to the real files/directories on the TSCC cluster
(the data are too large to store in git — the linked assemblies are 2–3 GB each).
They document *what* each stage reads and *where* it comes from.

## Data availability — download from the source databases

The data are not distributed with this repository. Please obtain them from the
public archives below.

| Data | Archive | Accession |
|---|---|---|
| HiFi long-read + Illumina short-read (this study) | NCBI BioProject | **PRJNA1393203** |
| 10x Single Cell Multiome | NCBI GEO | **GSE314504** |
| Hi-C | NCBI GEO | **GSE314503** |
| CENP-A and H3K27ac CUT&Tag | NCBI GEO | *processing; accessions to be assigned* |
| Comparative reference genomes (mouse, human, rat, dog, guinea pig, naked mole-rat) | NCBI RefSeq / GenBank | see `references/README.md` |

- NCBI BioProject: https://www.ncbi.nlm.nih.gov/bioproject/
- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/

While these accessions remain under embargo, access is granted through
reviewer links and GEO access tokens issued with the manuscript. Those
credentials are **deliberately not stored in this repository** — see the
manuscript's Data Availability section, or contact the corresponding author.

## Genome browser

An interactive genome browser for the OctDeg2.0 assembly and associated
features is available at:

https://renlab.sdsc.edu/share/jhc103/degu-genome-browser/

## To reproduce

1. Download the relevant BioProject / GEO data using the accessions above.
2. Either (a) replace these symlinks with your downloaded files, keeping the
   same names, or (b) edit `config.sh` (repo root) to point at your own copies.
3. `config.sh` is the single source of truth for all paths. Every script and
   notebook reads it — no script or notebook contains a hardcoded data path.

## Contents

| Directory | What | Key `config.sh` vars |
|---|---|---|
| `assemblies/` | de-novo assemblies (mito-assigned, masked, contam-filtered, haplotig-purged, hap2 scaffolds, verkko hap1, HapHiC build) | `MITO_ASSEMBLY`, `MASKED_ASSEMBLY`, `CHR_ASSIGNED_ASSEMBLY`, `FINAL_ASSEMBLY`, `VERKKO_HAP1_SCAFFOLDS_FA`, `HAPHIC_ASSEMBLY_FA` |
| `reads/` | HiFi (PacBio), Hi-C, RNA-seq (SRA), Illumina WGS | `HIFI_READS`, `HIC_*`, `MRNA_DIR`, `WGS_READS` |
| `annotations/` | de-novo (Braker) / Liftoff annotation GFFs | `LIFTOFF_REF_GFF`, `PARALOG_GFF` |
| `references/` | mouse GRCm39, OctDeg1, and comparative genomes | `MOUSE_GENOME_FA`, `HAPHIC_REF`, `HUMAN_GENOME_FA`, `RAT_GENOME_FA`, `DOG_GENOME_FA`, `GUINEAPIG_GENOME_FA`, `NMR_GENOME_FA` |

## Note

The symlinks target absolute paths on this cluster and will be broken after a
`git clone` elsewhere — that is expected. Git stores a symlink as its target
path only (mode `120000`), never the file contents, so **no sequence data is
committed by this directory**. Use the archived data plus `config.sh` to
reconstruct the inputs on your own system.
