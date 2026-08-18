# Inputs (stand-in)

This directory is a **stand-in** for the input data used by the pipeline. The
entries are **symbolic links** to the real files/directories on the TSCC cluster
(the data are too large to store in git). They document *what* each stage reads
and *where* it comes from.

## To reproduce

The real input data will be released through NCBI BioProjects — see the
manuscript's Data Availability section for accession numbers. To rerun:

1. Download the relevant BioProject files.
2. Either (a) replace these symlinks with the downloaded files, or (b) edit
   `config.sh` (repo root) to point at your own copies.
3. `config.sh` is the single source of truth for all paths; the symlinks here
   mirror its `*_ASM`, `*_READS`, `*_GFF`, and reference variables.

## Contents

| Directory | What | Key `config.sh` vars |
|---|---|---|
| `assemblies/` | de-novo assemblies (mito-assigned, masked, contam-filtered, haplotig-purged, hap2 scaffolds) | `MITO_ASSEMBLY`, `MASKED_ASSEMBLY`, `CHR_ASSIGNED_ASSEMBLY`, `FINAL_ASSEMBLY`, `REPEAT_MASKER_ASM` |
| `reads/` | HiFi (PacBio), Hi-C, RNA-seq (SRA), Illumina WGS | `HIFI_READS`, `HIC_*`, `MRNA_DIR`, `WGS_READS` |
| `annotations/` | de-novo (Braker) / Liftoff annotation GFFs | `LIFTOFF_REF_GFF`, `PARALOG_GFF` |
| `references/` | mouse GRCm39, OctDeg1 reference | `MOUSE_GENOME_FA`, `HAPHIC_REF` |

## Note

The symlinks target absolute paths on this cluster and will be broken after a
`git clone` elsewhere — that is expected. Use the BioProject data + `config.sh`
to reconstruct them on your own system.
