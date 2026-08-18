# 02-paralog-read-depth

BUSCO-calibrated HiFi read-depth evaluation of annotated paralog loci —
distinguishing genuine segmental duplications from uncollapsed haplotypes.

## Workflow (upstream → screen → notebook)

1. `build_paralog_families.py` — build paralog families (parent + `-lN`/`-dlN`/`-rlN` copies) from the annotation GFF → `paralog_families.{gff,tsv,_summary.tsv}`.
2. `read-depth-screen/` — map HiFi reads to the assembly and screen coverage across parent–paralog pairs.
   - `run_read_depth_screen.sh` (minimap2 + `read_depth_screen.py`)
3. `rna-seq-mapping/` — RNA-seq mapping + counting for paralog expression support.
   - `rna_seq_level.sh` (STAR index) → `run_star_alignment.sh` → `run_paralog_featureCounts.sh`
4. `notebook/category_exploration_v4_busco_used.ipynb` — BUSCO-calibrated read-depth evaluation (panels A–J), built from `notebook/build_category_exploration_nb.py`.

## Notes

- `read_depth_screen.py` and `build_category_exploration_nb.py` are argparse- or
  relative-path-driven; the `.sh` runners source `config.sh` (repo root).
- Data inputs mirror `config.sh`: `READ_DEPTH_ASSEMBLY`, `READ_DEPTH_HIFI_DIR`,
  `PARALOG_GFF`, `PARALOG_REPEAT_BED`, `PARALOG_GENOME_LENGTHS`, `STAR_INDEX`,
  `MRNA_DIR`.
