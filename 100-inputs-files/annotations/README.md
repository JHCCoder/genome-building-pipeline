# annotations/ — annotation GFFs (symlinks only)

The entries here are **symbolic links** to annotation GFF files on the TSCC
cluster. No annotation data is stored in this repository.

These are outputs of this study's annotation pipeline (stage 03: Liftoff
transfer + Braker3 de novo prediction, merged), not downloaded resources. To
regenerate them you need the assembly, the reference annotation, and the RNA-seq
evidence.

| Data | Archive | Accession |
|---|---|---|
| RNA-seq (Braker3 evidence) | NCBI BioProject | **PRJNA1393203** |
| 10x Single Cell Multiome | NCBI GEO | **GSE314504** |

- NCBI BioProject: https://www.ncbi.nlm.nih.gov/bioproject/
- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/

While these accessions remain under embargo, reviewer links and GEO access
tokens are issued with the manuscript. Those credentials are **not stored in
this repository**.

## Contents

| Link name | What | `config.sh` var |
|---|---|---|
| `denovoEnhanced_peaks2utr_sorted.gff3` | merged annotation with peaks2utr UTRs | `LIFTOFF_REF_GFF` |
| `denovoEnhanced_peaks2utr_sorted_perfect.gff` | earlier filtered variant (superseded) | — |

The current canonical annotation is `FINAL_ANNOTATION_GFF` in `config.sh`
(`..._parentReassigned_OptionB_20260820.gff3`, 26,273 genes). A gzipped copy of
the AGAT-normalized annotation is committed under
`../../100-outputs-genome-information/annotation/`.

An interactive genome browser showing these annotations is available at:

https://renlab.sdsc.edu/share/jhc103/degu-genome-browser/
