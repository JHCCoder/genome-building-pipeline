# List-building notebooks (add / rename / replace)

Three self-contained notebooks that build the three curated gene lists consumed
by `annotation-merging` (`merge_annotations.py`). Together they implement the
three enhancement scenarios used to complement a Liftoff-transferred annotation
with a de novo (Braker3) annotation.

Each notebook is available as a `.py` script (source of truth, runnable
headlessly) and a matching `.ipynb` notebook (regenerated with
`build_notebooks.py`).

## Run order

| # | Notebook | Scenario | Builds | Consumes |
|---|---|---|---|---|
| 1 | `02_rename_loc_genes` | 2 — rename | `gene_LOC_nameChange_unique.tsv` | — |
| 2 | `03_replace_loc_genes` | 3 — replace | `gene_LOC_replace_unique.tsv`, `gene_LOC_replace_unique_with_annotated_gene.tsv` | rename list (name-collision check) |
| 3 | `01_add_novel_genes` | 1 — add | `gene_gBraker_toAdd_unique.tsv` | rename list (Liftoff naming) |

The three lists are then passed to `merge_annotations.py` as
`--rename-map`, `--replace-list`, and `--add-list` (list mode).

## The three scenarios

1. **Add** — add a de novo gene absent from Liftoff when it (i) has no CDS
   overlap with any Liftoff gene, (ii) its CDS length is ≥ 60 % of its
   same-named mouse (mm39 / GRCm39) homolog, and (iii) it is expressed
   (BigWig signal > 0). Paralogs are suffixed `-dl1`, `-dl2`, ….
2. **Rename** — rename a `LOC` Liftoff gene to a de novo name when the two are
   the same locus (≥ 90 % same-strand CDS overlap) and the congruence is unique
   (1:1). Paralogs are suffixed `-l1`, `-l2`, ….
3. **Replace** — replace multiple `LOC` Liftoff genes with one de novo gene when
   that de novo gene is congruent with several `LOC` genes. Paralogs are
   suffixed `-rl1`, `-rl2`, ….

## Inputs

| Notebook | Liftoff GFF | De novo GFF | Other |
|---|---|---|---|
| add | `hifiasm-041425-scaffolded-chrAssigned-mito.gff` (genes), `gffl1.gff` (CDS) | `gffd2_sorted.gff` | mouse `GCF_000001635.27_GRCm39_genomic.gff`, `merged.bw` (BigWig), rename list |
| rename / replace | `hifiasm-041425-scaffolded-chrAssigned-mito_agatProcessed.gff` | `braker_3UTRincluded_agatProcessed.gff` | (rename list for replace) |

Paths are set in the `## 0. Configuration` cell of each notebook. The defaults
reproduce the *Octodon degus* annotation; edit them for your own data.

## Notes

- Output lists are written with a leading (unnamed) index column — the format
  expected by `merge_annotations.py` and the `01_add_novel_genes` notebook.
- The rename/replace notebooks require the input GFFs to be **sorted** (gene →
  mRNA → CDS/exon), which AGAT-processed GFFs are.
- `02` and `03` share the CDS-overlap step; each recomputes it (~minutes on a
  chromosome-scale genome). Keep them separate so each scenario is a standalone,
  auditable unit.
