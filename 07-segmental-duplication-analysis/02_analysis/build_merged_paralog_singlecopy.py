#!/usr/bin/env python3
"""Build the merged paralog + single-copy tables used by the expression figures.

`12_expression_violin_for_genes_of_interest.ipynb` and `12b_plot_figure5D_expression_violin.py` (Figure 5D) reads ONLY these two merged
tables -- not `paralog_families*` directly -- so they must be rebuilt whenever
the paralog tables or the annotation change, or the figure silently keeps the
previous gene names.

Inputs
  <PARA>/paralog_families_counts.tsv   featureCounts over the 3,272 genes of the
                                       multi-member families
  <PARA>/paralog_families.tsv          family metadata for the same genes
  ./interferon_singlecopy_counts.tsv   featureCounts over the 3 interferon genes
                                       that have no paralogs (Gbp1, Casp1,
                                       Syncrip). These are single-member
                                       families, so family-parent reassignment
                                       never renames them.

Outputs
  ./merged_paralog_singlecopy_counts.tsv
  ./merged_paralog_singlecopy_families.tsv

Previously this lived only as a shell heredoc recorded in code_log.txt; it is a
script now so the regeneration is reproducible.

Usage
    python build_merged_paralog_singlecopy.py [--para-dir DIR] [--force]
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

DEFAULT_PARA = ("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/"
                "figure/paralog-alignment-visualization")

# The three interferon genes with no paralogs; keys are featureCounts Geneids.
SINGLECOPY_NAMES = {
    "gene-Syncrip": "Syncrip",
    "gene-Gbp1": "Gbp1",
    "gene-Casp1": "Casp1",
}

COUNTS_OUT = "merged_paralog_singlecopy_counts.tsv"
FAMILIES_OUT = "merged_paralog_singlecopy_families.tsv"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--para-dir", default=DEFAULT_PARA)
    ap.add_argument("--singlecopy-counts", default="interferon_singlecopy_counts.tsv")
    ap.add_argument("--force", action="store_true",
                    help="Overwrite the outputs (they are regenerated in place).")
    args = ap.parse_args()

    for out in (COUNTS_OUT, FAMILIES_OUT):
        if os.path.exists(out) and not args.force:
            raise SystemExit(f"ERROR: {out} exists; pass --force to regenerate "
                             f"(back it up first).")

    pc = pd.read_csv(f"{args.para_dir}/paralog_families_counts.tsv", sep="\t", comment="#")
    sc = pd.read_csv(args.singlecopy_counts, sep="\t", comment="#")

    # featureCounts writes BAM column headers as full paths when invoked with
    # absolute paths; reduce to basenames so the two tables align.
    sc.columns = [c.split("/")[-1] if c.startswith("/") else c for c in sc.columns]
    missing = [c for c in pc.columns if c not in sc.columns]
    if missing:
        raise SystemExit(f"ERROR: single-copy table is missing columns present in "
                         f"the paralog table: {missing[:5]}")
    sc = sc[pc.columns]

    merged_counts = pd.concat([pc, sc], ignore_index=True)
    merged_counts.to_csv(COUNTS_OUT, sep="\t", index=False)

    pf = pd.read_csv(f"{args.para_dir}/paralog_families.tsv", sep="\t")
    unknown = [g for g in sc["Geneid"] if g not in SINGLECOPY_NAMES]
    if unknown:
        raise SystemExit(f"ERROR: unrecognised single-copy Geneid(s): {unknown}. "
                         f"Update SINGLECOPY_NAMES.")

    rows = []
    for _, r in sc.iterrows():
        gid = r["Geneid"]
        rows.append(dict(family=SINGLECOPY_NAMES[gid], gene_name=SINGLECOPY_NAMES[gid],
                         gene_id=gid, gene_type="singlecopy", chrom=r["Chr"],
                         start=int(r["Start"]), end=int(r["End"]), strand=r["Strand"],
                         length=int(r["End"]) - int(r["Start"]) + 1,
                         paralog_type=None, copy_num=None))
    merged_families = pd.concat([pf, pd.DataFrame(rows)], ignore_index=True)
    merged_families.to_csv(FAMILIES_OUT, sep="\t", index=False)

    print(f"  {COUNTS_OUT}:   {len(merged_counts):,} rows "
          f"(paralog {len(pc):,} + single-copy {len(sc)})", file=sys.stderr)
    print(f"  {FAMILIES_OUT}: {len(merged_families):,} rows "
          f"(paralog {len(pf):,} + single-copy {len(rows)})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
