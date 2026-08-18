# 100-outputs-genome-information

Actual output / result files for the genome assembly, kept separate from the
pipeline code (stages `01–07`). Organized by analysis tool. Only small,
text-based reports are stored here; large files (assemblies, read alignments,
raw reads) are deposited in a data repository (NCBI / Zenodo) and referenced by
accession rather than committed to git.

## Contents

| Directory | Analysis | Tool |
|---|---|---|
| `Inspector/` | assembly structural-error detection | Inspector v1.3.1 |

### `Inspector/`

Three runs, each with the report files `summary_statistics`,
`structural_error*.bed`, `small_scale_error*.bed`, `Inspector.log`, and
`contig_length_info`:

- `hifiasm-041425-scaffold/` — the scaffolded, chromosome-assigned, mito-added
  assembly (the final OctDeg2.0).
- `hifiasm-041425-contig/` — the pre-scaffolding contig assembly.
- `hifiasm-041425-contig-reference/` — the contig assembly, reference-based mode
  (adds `structural_errors_ref.bed` / `small_scale_error_ref.bed` against OctDeg1).

The full Inspector outputs (read-to-contig BAMs, `valid_contig.fa` assemblies,
~148 GB) are **not** stored here — only the ~25 MB of report files.
