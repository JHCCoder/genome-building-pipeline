# 04-paralog-expression

Paralog RNA-seq expression analysis: Salmon transcript quantification and the
paralog-expression heatmaps / support analysis.

## Contents

- `salmon-mapping/` — Salmon transcriptome quantification:
  - `00_building_index.sh` (index) → `01_mapping_reads.sh` (quant) → `02_merge_results.sh` (quantmerge → TPM matrix)
- `make_paralog_sample_heatmap_R_used.ipynb` (+ `make_paralog_sample_heatmap.py`/`.R`) — expression heatmap across 29 tissues (ComplexHeatmap).
- `paralog_support_analysis.ipynb` — parent–paralog RNA-seq support classification.

## Inputs

- `paralog_families_counts.tsv` — featureCounts / Salmon counts (from `03-paralog-read-depth/` or `salmon-mapping/`).
- `paralog_families.tsv` — family metadata (from `build_paralog_families.py` in `03-paralog-read-depth/`).

## Notes

- The notebooks and `.py`/`.R` read their inputs as relative paths (working
  directory); the Salmon `.sh` runners source `config.sh` (repo root) for the
  assembly/GFF/read paths and conda env.
