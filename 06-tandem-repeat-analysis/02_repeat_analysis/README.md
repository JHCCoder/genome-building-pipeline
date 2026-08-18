# 05-tandem-repeat-analysis — analysis

Analysis/visualization of the tandem-repeat results (Supplemental Figure S15).

## Notebooks (S15 panels)

| Panel | Notebook | What it shows |
|---|---|---|
| A | `check_tandem_repeats_degu.ipynb` | t-SNE of longest tandem-repeat consensus (degu, 3 period classes, 83 sequences) |
| B | `top_tandem_repeat_multiple_species.ipynb` | clustalo tree (10 longest repeats × 7 species: degu, guinea pig, mouse, NMR, rat, dog, human) |
| C | `check_tandem_repeats_human.ipynb` | human TRF scatter |
| D | `check_tandem_repeats_mouse.ipynb` | mouse TRF scatter |
| E | `check_tandem_repeats_NMR.ipynb` | naked mole-rat TRF scatter |
| — | `check_tandem_repeats_dog/guineaPig/rat.ipynb`, `check_tandem_repeats_degu_verkko.ipynb` | remaining species for the tree / alternate assembly |

## Classification scripts (`classification/`)

- `classify_trf_hits_caseinsensitive.py` — classify TRF hits into period classes.
- `classifier_to_igv_bed.py`, `171peak_tsv_to_igv_bed.py`, `ucsc_mm39_gap_to_igv_bed.py` — convert classified hits to IGV BED tracks.

## Supplemental Figure S19

- `top10_HOR_lengths_per_chromosome.ipynb` — per-chromosome count distribution
  of the top-10 HOR lengths (centroAnno `chr*_horDecomposedResult.csv`), coloured
  by average identity.
- `top10_monomer_lengths_per_chromosome.ipynb` — same, for the top-10 monomer
  lengths (`chr*_decomposedResult.csv`).

## Notes

- These read the TRF `.dat` / classified outputs (produced by the `01-trf/`
  processing batch) from the working directory — see `../processing/`.
- The two degu notebooks are large (embedded plot outputs). `check_tandem_repeats_degu.ipynb`
  here is the canonical `code/command-line-script/...` version (per FIGURE_CODE_MAP).
