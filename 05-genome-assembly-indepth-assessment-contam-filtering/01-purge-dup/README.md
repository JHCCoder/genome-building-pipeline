# 01-purge-dup

Purge haplotypic duplicates from the scaffolded assembly with `purge_dups`,
plus the analysis behind **Supplemental Figure S13** ("One additional round of
purged duplication does not substantially affect genome annotation").

## Command pipeline (step1–4)

| Step | Script | What it does |
|---|---|---|
| 1 | `step1_align.sh` | minimap2 read alignment + `pbcstat`/`calcuts` + self-alignment |
| 2 | `step2_purge_haplotids_overlaps.sh` | `purge_dups` → `dups.bed` |
| 3 | `step3_purge_from_assembly.sh` | `get_seqs` (remove flagged dups) |
| 4 | `step4_round2.sh` | second round (hap.fa + haplotig assembly) |

## S13 figure panels

| Panel | Content | Code |
|---|---|---|
| A | purge_dup results table (HAPLOTIG / OVLP / HIGHCOV / JUNK / REPEAT) | `dups.bed` |
| B | Genome-wide distribution of purge_dup + FCS-GX | `purge-dup-distribution/purge_dup_distribution.ipynb` |
| C | De-novo genes/paralogs per purge_dup category | `paralog-gene-overlap/haplotig_gene_overlap.ipynb` (+ `recompute_haplotig_overlap.py`) |
| D | 12.09 Mb haplotig distribution per chromosome | `purge-dup-distribution/haplotig_scaffold_distribution.ipynb` |
| E | REPEAT vs BISER segdup overlap vs 10,000-permutation null | `biser-enrichment/` |
| F | BISER segdup vs purge_dups category overlap | `biser-enrichment/` |

Only panel E involves segmental duplication (BISER); the rest are purge_dup-only.

## Notes

- `biser-enrichment/` is a self-contained project (`config/` + `scripts/` +
  `notebook/`); its Python scripts hardcode `PROJECT_ROOT` — edit it to match
  `config.sh` (repo root).
- All notebooks read their inputs from a clearly-marked configuration cell at
  the top (mirrors `config.sh`).
