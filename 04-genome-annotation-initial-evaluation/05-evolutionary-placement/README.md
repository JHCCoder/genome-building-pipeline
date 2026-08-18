# 05-evolutionary-placement

Phylogenetic placement of the degu among related rodent species, reconstructed
from the methods section (the original analysis code was not retained). The
pipeline finds conserved single-copy orthologs with BUSCO, aligns and trims each
ortholog with MAFFT + trimAl, concatenates them into an amino-acid supermatrix,
and infers a maximum-likelihood tree with RAxML.

## Pipeline steps (in order)

| # | Step | Script(s) |
|---|------|-----------|
| 1 | BUSCO (genome mode, `euarchontoglires_odb12`, `--metaeuk`, `-c 16`) per species | `01_busco_genome_mode.sh` |
| 2 | Retain complete single-copy genes present in *every* genome; write per-gene multi-species FASTAs | `02_collect_shared_singlecopy_orthologs.sh` (+ `.py`) |
| 3 | Align (MAFFT) + trim (trimAl) each gene independently | `03_align_trim_concat.sh` |
| 4 | Concatenate trimmed alignments into one PHYLIP supermatrix | `04_concat_supermatrix.py` (invoked by step 3) |
| 5 | RAxML maximum-likelihood tree (JTT+Γ, `-f a -N 100 -T 16`) | `05_raxml_phylogeny.sh` |

The methods-equivalent BUSCO command reproduced by step 1 is:

```
busco -i <genome.fna> -m genome -l euarchontoglires_odb12 -c 16 --metaeuk
```

## Inputs

- `species_list.tsv` — one `<genome.fna>	<tree_label>` per line. The default
  set is the degu assembly plus the related rodents and dog/human outgroups
  listed in `data/genome-related-species/genome_list6.txt`. Edit this file to
  change which genomes are included; keep labels ≤ 10 chars, unique, and
  whitespace-free (RAxML relaxed-PHYLIP name limit).

## Outputs (under `$EVO_OUT_DIR`)

| Directory | Contents |
|-----------|----------|
| `01_busco/` | one BUSCO run per species (`<label>/run_<lineage>/`) |
| `02_single_copy_orthologs/` | `<gene>.faa` — one multi-species FASTA per shared gene |
| `03_alignments/` | `<gene>.aln.faa` (MAFFT) and `<gene>.trim.faa` (trimAl) |
| `04_supermatrix/` | `supermatrix.phy` — concatenated relaxed-PHYLIP supermatrix |
| `05_raxml/` | `RAxML_bipartitions.degu_tree` — final tree with bootstrap support |

## Running

Run the numbered scripts in order. Steps 1, 3 and 5 are Slurm batch jobs; step 2
is a lightweight script you run directly after step 1 finishes:

```bash
# Step 1 — BUSCO on every species genome
sbatch 01_busco_genome_mode.sh

# Step 2 — collect shared single-copy orthologs (after step 1 finishes)
bash 02_collect_shared_singlecopy_orthologs.sh

# Step 3 — MAFFT + trimAl per gene, then concatenate (step 4)
sbatch 03_align_trim_concat.sh

# Step 5 — RAxML tree
sbatch 05_raxml_phylogeny.sh
```

## Configuration

Paths, conda environments, and the BUSCO lineage live in `config.sh` at the
repository root (section "04-genome-annotation-initial-evaluation —
evolutionary-placement"). Each shell script sources it automatically; edit
`config.sh` (not the scripts) for your environment.

## Requirements

- `toolshed-busco` conda env — BUSCO + MetaEuk + Augustus (used by step 1).
- `evolutionary-tree` conda env — MAFFT 7.526 and RAxML 8.2.12, plus trimAl 1.5.
  trimAl is not always preinstalled; if `trimal` is missing after
  `conda activate evolutionary-tree`, install it with
  `conda install -c bioconda trimal`.
- Python 3 (steps 2 and 4 use only the standard library).

## Notes

- The methods cite BUSCO 5.8.2; the `toolshed-busco` env ships BUSCO 5.7.1.
  Outputs are equivalent for this use; re-point `ENV_BUSCO` if 5.8.2 is required.
- `euarchontoglires_odb12` is auto-downloaded by BUSCO on the first run if it is
  not already cached locally.
- RAxML's `-p`/`-x` random seeds are not stated in the methods; fixed seeds
  (12345/12345) are used for reproducibility and can be changed at the top of
  `05_raxml_phylogeny.sh`.
- trimAl's `-automated1` and MAFFT's `--auto` are sensible defaults for
  single-copy-ortholog phylogenomics; both are easy to override in
  `03_align_trim_concat.sh`.
- Scripts are Slurm batch scripts written for TSCC at UCSD; adjust the `#SBATCH`
  scheduler lines for your own cluster.
