# 07-HiC-contact-map-visualization

Hi-C contact-map visualization pipeline (HiCExplorer + cooltools). Generates Hi-C
matrices from aligned Hi-C reads, plots whole-genome / per-chromosome contact maps,
and calls A/B compartments.

Pulled from `figure/hic-plot/` and generalized for the shared repo (paths and conda
environments live in `config.sh`; scheduler email generalized).

## Pipeline steps (in order)

| # | Step | Script |
|---|------|--------|
| 1 | Align HiFi Hi-C reads (bwa-mem2 → BAM) | `01_align_for_hicexplorer.sh` |
| 2 | Align deep-file Hi-C reads | `02_align_for_hicexplorer_deepFile.sh` |
| 3 | Short-read (Illumina) Hi-C: align → matrix → merge | `03_hicexplorer_shortread_pipeline.sh` |
| 4 | Build Hi-C matrix (`hicBuildMatrix`, 10 kb) | `04_build_matrix.sh` |
| 5 | Per-chromosome matrix + A/B compartments (`hicPCA`) | `05_build_matrix_per_chromosome.sh` |
| 6 | Merge matrix bins (`hicMergeMatrixBins`) | `06_merge_matrix_bins.sh` |
| 7 | Plot contact map (`hicPlotMatrix`) | `07_plot_matrix.sh` |
| 8 | Split BAMs by chromosome (helper) | `08_split_chromosome.sh` |
| 9 | Submit the split-by-chromosome array (helper) | `09_split_chromosome_submission.sh` |
| 10 | A/B compartment eigenvector (cooltools) | `10_ab_compartment_calculation.ipynb` |

## Inputs (paths set in `config.sh`)

| Variable | Data |
|----------|------|
| `HIC_GENOME` | `…/assembly_final.sorted.headerRenamed.chrAssigned.fasta` (alignment + restriction sites) |
| `HIC_RESTRICTION_BED` | `figure/hic-plot/041425_assembly_GATC.bed` (GATC sites) |
| `HIC_HIFI_READ_DIR` | trimmed HiFi Hi-C fastq (`haphic/input-hic-read`) |
| `HIC_SHORTREAD_DIR` / `$DATA_DIR/sequencing-reads-HiC` | Illumina / raw Hi-C fastq |
| cooler `.mcool` | `output/matrix/403_ear_deep/403_ear_deep_octDeg2.mcool` (notebook) |

The restriction-site BED is regenerated with:

```
hicFindRestSite --fasta "$HIC_GENOME" --searchPattern GATC -o "$HIC_RESTRICTION_BED"
```

## Conda environments (in `config.sh`)

- `ENV_BULK_HIC` (`bulk-HiC-processing`) — bwa-mem2, trim_galore, samtools
- `ENV_HICEXPLORER` (`toolshed-HiCExplorer`) — hicBuildMatrix / hicPlotMatrix / hicPCA / hicMerge*
- `ENV_SAMTOOLS` (`toolshed-samtools`) — samtools + sambamba
- Notebook needs Python with `cooler`, `cooltools`, `bioframe`, `pandas`, `numpy`, `matplotlib`.

## Outputs (under `figure/hic-plot/`)

- `hicMatrix/`, `hicMatrix-longRead/` — `.h5` matrices
- `plots/` — `hicPlotMatrix` PNGs
- `compartments/` — per-chromosome A/B compartment BEDs
- `eigenvector_track.tsv` / `eigenvector_track_dropNA.bed` — A/B eigenvector track

## Notes

- Scripts are Slurm batch scripts written for TSCC at UCSD; adjust the `#SBATCH`
  scheduler lines for your own cluster.
- Step 5 (`05_build_matrix_per_chromosome.sh`) currently hardcodes `chr1 chr2` in the
  `CHROMOSOMES` array; edit it (or use `chromosome_list.txt`) for the full genome.
