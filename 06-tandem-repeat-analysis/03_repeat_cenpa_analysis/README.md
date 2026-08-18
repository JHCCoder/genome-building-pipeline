# 03_repeat_cenpa_analysis

CENP-A CUT&Tag and centromeric repeat analyses (Figure 4, Figure 5, S17–S19).

## Subdirectories

| Directory | Figures | Content |
|---|---|---|
| `cenpa-cuttag-enrichment/` | Fig 4, Fig 5B–C | CENP-A CUT&Tag enrichment vs TRF period bins. `figure_notebook.ipynb`, `prep_hicat_karyotype.py`, numbered pipeline (`1_validate_and_prepare.sh` … `5_analysis.R`), karyotype `*.R`, and sub-analyses (`period-enrichment/`, `349-bp/`, `195bp-repeatmasker-overlap/`). |
| `genome-wide-visualization/` | Fig 5A | Genome-wide telomere/satellite/HOR/CENP-A tracks (`genomewide_telomere_track.R`, `karyotype_centroAnno_telomere.R`, `scan_telomeres.py`, `analyze_telomere_overlap.py`, `trf_telomere_array.sh`). |
| `circos/` | Fig 5A (related) | CENP-A circos (`pyCircularize_CENPA_195_389_349.ipynb`). |
| `cenpa-repeat-chromosome/` | S17 | TRF-array occupancy around CENP-A peaks (`repeat_occupancy_genomewide_plot.ipynb`, `repeat_occupancy_around_cenpa_core.ipynb` + numbered `scripts/` pipeline). |
| `centromeric-HOR-metrics/` | S18, S19 | centroAnno monomer / HOR length decomposition (`monoDecomposeResult.ipynb`, `horDecomposeResult.ipynb`, + 3 supplementary notebooks). |
| `stained-glass/` | Fig 5D | StainedGlass runner (`submit_stainedGlass.sh`) — third-party Snakemake tool, not vendored. |

## Notes

- **StainedGlass** (Fig 5D) is a third-party Snakemake tool (`toolshed/StainedGlass/`);
  only the degu-specific runner `stained-glass/submit_stainedGlass.sh` is included,
  not the tool itself.
- Exploratory `*_tmp.py` scripts and `*_backup*.ipynb` notebooks were excluded.
- Large notebooks had their embedded outputs stripped; re-run them to regenerate figures.
