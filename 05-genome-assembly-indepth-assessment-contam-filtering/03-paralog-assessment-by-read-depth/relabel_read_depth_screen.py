#!/usr/bin/env python3
"""
Relabel the read-depth-screen outputs to the parent-reassigned annotation.

The family-parent reassignment (2026-08-20) is a PURE RENAME: every gene keeps
its (chrom, start, end) but 708 genes across 268 families receive a new Name/ID
and a possibly-swapped parent/copy role.  Read-depth *values* are keyed to
coordinates and are therefore unchanged — only the gene-name / gene-type labels
change.  The locus table is relabeled in place; the pair table is REBUILT (not
just relabeled) because `build_pair_table` forms parent–copy pairs from the
`gene_type` column, which the reassignment changes.

Inputs (under PROJ_ROOT):
  family_parent_rename_map_20260820.tsv            old/new name+id, role (all 3,272 family members)
  figure/paralog-alignment-visualization/paralog_families.tsv   NEW authoritative gene_type/paralog_type/copy_num
  figure/…/read-depth-screen-v4-busco-calibrated/locus_classification_v4_busco.tsv  (relabel)
  figure/…/read-depth-screen-v4-busco-calibrated/pair_integration_v4_busco.tsv       (rebuild)
  figure/paralog-alignment-visualization/paralog_genomic_summary.tsv  (relabel; notebook input #3)
  figure/paralog-alignment-visualization/paralog_support_summary.tsv  (relabel; feeds pair-table expression columns)
  figure/…/read-depth-screen-v4-busco-calibrated/selected_busco_thresholds.tsv       (reconstruct ThresholdSet)
  figure/…/read-depth-screen-v4-busco-calibrated/calibrate_and_classify_v4.py         (import build_pair_table + make_threshold_set)

Outputs overwrite three files in place (after a timestamped backup):
  - locus_classification_v4_busco.tsv   (relabeled)
  - pair_integration_v4_busco.tsv       (rebuilt)
  - paralog_genomic_summary.tsv         (relabeled; notebook input #3)

Run with --dry-run first to verify; pass --write to commit.

Note: `paralog_support_summary.tsv` (expression) is used only IN MEMORY to carry
the pair table's expression columns; it is NOT written back.  It is doubly stale
(2,248 rows from the pre-regeneration annotation; 223 names have no match in the
current rename map) and belongs to the Supp Fig S9 (expression) regeneration.
The pair table's expression columns (expression_category / allele_risk /
segdup_support) are carried best-effort and are NOT consumed by the S8 notebook,
which uses only pair_interpretation / locus_coverage_category / paralog_type /
gene_type / gene_name / family.
"""

import argparse
import importlib.util
import os
import shutil
import sys

PROJ_ROOT = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"

RENAME_MAP = os.path.join(
    PROJ_ROOT,
    "code/command-line-script/annotation-merging/output/"
    "family_parent_rename_map_OptionB_from_expr_20260820.tsv",
)
FIG = os.path.join(PROJ_ROOT, "figure/paralog-alignment-visualization")
SCREEN = os.path.join(FIG, "read-depth-screen-v4-busco-calibrated")

PARALOG_FAMILIES = os.path.join(FIG, "paralog_families.tsv")
LOCUS = os.path.join(SCREEN, "locus_classification_v4_busco.tsv")
PAIR = os.path.join(SCREEN, "pair_integration_v4_busco.tsv")
GENOMIC = os.path.join(FIG, "paralog_genomic_summary.tsv")
SUPPORT = os.path.join(FIG, "paralog_support_summary.tsv")
THRESHOLDS = os.path.join(SCREEN, "selected_busco_thresholds.tsv")
CALIBRATE = os.path.join(SCREEN, "calibrate_and_classify_v4.py")


def load_calibrate():
    spec = importlib.util.spec_from_file_location("calib", CALIBRATE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_threshold_set(calib):
    """Reconstruct the V4_BUSCO ThresholdSet from selected_busco_thresholds.tsv."""
    import pandas as pd

    t = pd.read_csv(THRESHOLDS, sep="\t").iloc[0].to_dict()
    # make_threshold_set expects exactly these keys; the tsv has extra columns
    # which are harmless (dict.get).
    return calib.make_threshold_set("V4_BUSCO", t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true",
                    help="overwrite outputs in place (default: dry-run only)")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    calib = load_calibrate()
    ts = build_threshold_set(calib)

    # --- load inputs -------------------------------------------------------
    rename = pd.read_csv(RENAME_MAP, sep="\t")
    fams = pd.read_csv(PARALOG_FAMILIES, sep="\t")
    locus = pd.read_csv(LOCUS, sep="\t")
    genomic = pd.read_csv(GENOMIC, sep="\t")
    support = pd.read_csv(SUPPORT, sep="\t")

    # --- build id/name/type mappings --------------------------------------
    # rename map: old_gene_id -> new_gene_id, new_gene_name  (covers all 3,272)
    old_id_to_new = dict(zip(rename["old_gene_id"], rename["new_gene_id"]))
    old_name_to_new = dict(zip(rename["old_name"], rename["new_name"]))

    # new paralog families: gene_id -> (gene_name, gene_type, paralog_type, copy_num)
    fam_by_id = {
        r["gene_id"]: (r["gene_name"], r["gene_type"],
                       r["paralog_type"], r["copy_num"])
        for _, r in fams.iterrows()
    }

    # --- relabel locus table ----------------------------------------------
    n_locus = len(locus)
    unmapped_locus = []
    for idx, row in locus.iterrows():
        oid = row["gene_id"]
        nid = old_id_to_new.get(oid)
        if nid is None or nid not in fam_by_id:
            unmapped_locus.append(oid)
            continue
        name, gtype, ptype, cnum = fam_by_id[nid]
        locus.at[idx, "gene_id"] = nid
        locus.at[idx, "gene_name"] = name
        locus.at[idx, "gene_type"] = gtype
        # paralog_type / copy_num: NaN for parents (matches original schema)
        locus.at[idx, "paralog_type"] = ptype if ptype else np.nan
        locus.at[idx, "copy_num"] = cnum if not pd.isna(cnum) else np.nan

    # --- relabel sidecar tables (paralog_name column) ----------------------
    def relabel_name_col(df, col):
        n_unmapped = 0
        for idx, row in df.iterrows():
            old = row[col]
            new = old_name_to_new.get(old)
            if new is None:
                n_unmapped += 1
            else:
                df.at[idx, col] = new
        return n_unmapped

    n_gen_unmapped = relabel_name_col(genomic, "paralog_name")
    n_sup_unmapped = relabel_name_col(support, "paralog_name")

    # --- rebuild pair table ------------------------------------------------
    pairs = calib.build_pair_table(locus, ts, support, genomic)

    # --- verify ------------------------------------------------------------
    print("=" * 70)
    print("VERIFICATION")
    print("=" * 70)
    print(f"locus rows          : {n_locus}  (expect 3,272)")
    print(f"locus unmapped ids  : {len(unmapped_locus)}  (expect 0)")
    if unmapped_locus:
        print("  first unmapped:", unmapped_locus[:10])
    print(f"locus gene_type     : {dict(locus['gene_type'].value_counts())}")
    print(f"genomic paralog_name unmapped : {n_gen_unmapped}  (expect 0)")
    print(f"support paralog_name unmapped : {n_sup_unmapped}  (expect 0)")

    n_pairs = len(pairs)
    n_copies = int((locus["gene_type"] != "parent").sum())
    print(f"pair rows           : {n_pairs}  (expect {n_copies} = one per copy)")
    print(f"pair interpretation :")
    for k, v in pairs["pair_interpretation"].value_counts().items():
        print(f"    {k}: {v}")

    # column parity with the original pair table
    orig_pair = pd.read_csv(PAIR, sep="\t")
    print(f"pair columns match original: "
          f"{list(pairs.columns) == list(orig_pair.columns)}")
    if list(pairs.columns) != list(orig_pair.columns):
        print("  new columns:", list(pairs.columns))
        print("  old columns:", list(orig_pair.columns))

    # spot-check a reassigned family (Gbp5): parent must now be the AUGUSTUS model
    gbp = pairs[pairs["family"] == "Gbp5"]
    print(f"\nGbp5 pairs (new parent/copy):")
    for _, r in gbp.iterrows():
        print(f"    {r['parent_gene']:>10}  <-  {r['paralog_gene']}  "
              f"[{r['pair_interpretation']}]")

    if not args.write:
        print("\nDry-run only. Re-run with --write to commit.")
        return 0

    # --- backup + write ----------------------------------------------------
    # support_summary is intentionally NOT written (see docstring).
    for path in (LOCUS, PAIR, GENOMIC):
        shutil.copy2(path, path + ".pre-relabel_20260821.bak")

    locus.to_csv(LOCUS, sep="\t", index=False)
    pairs.to_csv(PAIR, sep="\t", index=False)
    genomic.to_csv(GENOMIC, sep="\t", index=False)

    print("\nWrote (backups at *.pre-relabel_20260821.bak):")
    for path in (LOCUS, PAIR, GENOMIC):
        print("  " + path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
