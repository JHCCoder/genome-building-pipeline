# 06-mashTree-syntenic-evaluation

Synteny plot of the degu assembly (OctDeg2.0) against other rodent species
(**Supplemental Figure S4B**). The species are ordered along the plot's y-axis by a
**mashtree** built from the Jaccard distance between 21-mer sets.

Two independent sources of the species order are computed:

- **mashtree** (the one used in the figure) — a dendrogram built directly from
  21-mer Jaccard distances, in jackknife and bootstrap variants
  (`01_mashtree.sh`).
- **mash triangle → FastME** — an explicit pairwise Jaccard distance matrix
  turned into a BIONJ tree (`02_mash_distance_tree.sh`), included as an
  alternative distance tree.

The synteny itself is computed with ntSynt and drawn with ntSynt-viz.

## Order of operation

| Step | File | What it does |
|---|---|---|
| 1 | `01_mashtree.sh` | build jackknife + bootstrap mashtrees → `mashtree.{jackknife,bootstrap}.renamed.dnd` |
| 2 | `02_mash_distance_tree.sh` | (alternative) `mash sketch` → `mash triangle` → FastME BIONJ distance tree |
| 3 | `03_ntSynt.sh` | pairwise synteny blocks → `rodent.synteny_blocks.tsv` + per-genome `.fai` |
| 4 | `04_ntsynt_viz.sh` | ribbon plot with species ordered by the mashtree (`--tree`) |

Steps 1 and 2 are independent; step 4 depends on step 3 (blocks + `.fai`) and
step 1 (the mashtree).

## Inputs

- **`fasta_list.txt`** — one genome FASTA path per line: the degu assembly
  (`assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta`) plus 10 other
  genomes (8 rodents + human + dog outgroups) from `data/genome-related-species/`.
- **`name-conversions.tsv`** — maps raw FASTA filenames to display species names
  used on the plot. Two names carry typos from the source figure and are kept
  verbatim so the plot matches the manuscript: `Egytian_spiny_rat` (→ *Acomys*,
  "Egyptian") and `PPocket_mouse` (→ *Perognathus*, "Pacific pocket mouse").

Both files use absolute paths / fixed names and must be edited for a new
environment, or regenerated from the released assemblies.

## Configuration

Paths, tools, and environments live in `config.sh` (repository root). This stage
adds: `ENV_PHYLO`, `MASH_BIN`, `MASHTREE_DIR`, `NTSYNT_BIN`, `NTSYNT_VIZ_BIN`,
`SYNTENY_OUT_DIR`, `SYNTENY_FASTA_LIST`, `SYNTENY_NAME_CONVERSIONS`,
`SYNTENY_TARGET_GENOME`; it reuses `ENV_ASSESSMENT`.

## Requirements

- **mash** (`MASH_BIN`, standalone) + **mashtree** Perl scripts (`MASHTREE_DIR/bin`).
  mashtree needs Perl **BioPerl**, provided by the `evolutionary-tree` conda env
  (`ENV_PHYLO`), which also supplies **FastME**.
- **ntSynt** and **ntSynt-viz** (`NTSYNT_BIN`, `NTSYNT_VIZ_BIN`); ntSynt-viz runs
  a snakemake backend, provided by the `genome-assembly-assessment` env
  (`ENV_ASSESSMENT`).
- `rename_tree.py` / `convert_to_square.py` need Python 3 (`convert_to_square.py`
  additionally needs NumPy).

## Notes

- The scripts are Slurm batch scripts for TSCC/UCSD; adjust the `#SBATCH` lines
  for your cluster.
- The `.fai` files consumed by `04_ntsynt_viz.sh` are the ones written by ntSynt
  in step 3 (same genome names as the blocks TSV).
- Generated outputs (`.msh`, `.tab`, `.phy`, `.nwk`, `.dnd`, `.fai`,
  `.synteny_blocks.tsv`, `.png`) are large/intermediate and are not committed.
