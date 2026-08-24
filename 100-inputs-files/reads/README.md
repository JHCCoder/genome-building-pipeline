# reads/ — sequencing data (symlinks only)

The entries here are **symbolic links** to read directories on the TSCC cluster.
No sequencing data is stored in this repository.

## Download from the source archives

| Link | Data | Archive | Accession |
|---|---|---|---|
| `hifi/` | PacBio HiFi long reads (this study) | NCBI BioProject | **PRJNA1393203** |
| `illumina/` | Illumina short-read WGS (this study) | NCBI BioProject | **PRJNA1393203** |
| `hic/` | Hi-C | NCBI GEO | **GSE314503** |
| `rnaseq/` | RNA-seq (SRA deposits) | NCBI BioProject | **PRJNA1393203** |

- NCBI BioProject: https://www.ncbi.nlm.nih.gov/bioproject/
- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/

Related data from this study, not linked here: 10x Single Cell Multiome
(GEO **GSE314504**); CENP-A and H3K27ac CUT&Tag (GEO, currently processing —
accessions to be assigned).

While these accessions remain under embargo, reviewer links and GEO access
tokens are issued with the manuscript. Those credentials are **not stored in
this repository**.

## To reproduce

Download the reads, then either replace these symlinks with your own copies or
point the matching `config.sh` variables (`HIFI_READS`, `WGS_READS`, `HIC_*`,
`MRNA_DIR`) at them. `config.sh` is the single source of truth for all paths.
