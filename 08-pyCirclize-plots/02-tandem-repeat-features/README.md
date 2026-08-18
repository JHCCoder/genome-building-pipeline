# 02-tandem-repeat-features

pyCirclize circos plots of the degu tandem-repeat features. Pulled from
`figure/circos-plot/repeat-features/` and generalized for the shared repo
(paths use the `PROJ_ROOT` variable in the first code cell; tutorial and
exploratory cells removed; outputs cleared).

## Notebooks

| Notebook | Figure | Contents |
|----------|--------|----------|
| `01_repeat_features_circos.ipynb` | `degus_genome_circos_repeat_overview.png/.svg` | 349-bp centromeric repeats, HORs, tandem-repeat density, interspersed-repeat density |
| `02_cenpa_195_389_349_circos.ipynb` | `degus_genome_circos_CENPA_195_389_349.png/.svg` | CENP-A CUT&Tag (log10 bins) + 195/389-bp and 349-bp tandem repeats |

## Inputs (paths set via `PROJ_ROOT` in the first code cell)

| Data | Path |
|------|------|
| Chromosome sizes (`.fai`) | `data/…/assembly_final.sorted.headerRenamed.chrAssigned.hardMasked.mito.fasta.fai` |
| 349-bp centromeric repeats | `code/command-line-script/genome-annotation/trf-tandem-repeat/349peak_repeat_1millionBpMin.tsv` |
| 195/389-bp tandem repeats | `code/command-line-script/genome-annotation/trf-tandem-repeat/{195,389}peak_repeat.tsv` |
| HORs | `output/outputs-from-centraAnno/hifiasm-0414/cautils-chrOnly/HORs.bed` |
| Interspersed-repeat density | `figure/circos-plot/repeat-features/assembly_final.sorted.headerRenamed.repeatDensity.1mb_500kStep_rollingWindow.tsv` |
| Tandem-repeat density | `figure/circos-plot/repeat-features/assembly_final.sorted.headerRenamed.tanRepeatDensity.merged.tsv` |
| Scaffold junctions | `figure/circos-plot/feature-overview/agp_final_contig2scaffold.bed` |
| CENP-A bigWig | `figure/cenpa-cuttag-centromere/bw_files/XG_150.all.bw` |

## Requirements

Python with `pycirclize`, `biopython`, `pyBigWig` (notebook 2 only), `pandas`,
`numpy`, `matplotlib`.
