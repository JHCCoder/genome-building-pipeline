#!/usr/bin/env python3
"""Merge a Liftoff-transferred annotation with a de-novo (e.g. Braker) annotation.

The merge follows three enhancement rules (see README.md):

  1. name-change : optionally rename Liftoff genes from a mapping file.
  2. replace     : replace a Liftoff gene with a de-novo gene only when they
                   are the same gene (concordant name) or the Liftoff gene is
                   an unannotated LOC; keep differently-named Liftoff genes
                   (protects reference-derived annotations on paralog loci).
  3. add         : add functionally-named de-novo genes whose CDS does not
                   overlap any retained Liftoff gene.

Example
-------
    python merge_annotations.py \
        --liftOff   liftoff.gff \
        --denovo    braker.functional.gff \
        --outdir    output/ \
        --prefix    merged

Only the Python standard library is required (pyBigWig is optional and only
needed when an expression track is supplied).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, List, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.gffio import parse_gff, write_gff
from src.merge import MergeConfig, merge_annotations, merge_with_lists


def _strip_index(row: List[str]) -> List[str]:
    """Drop a leading pandas-style index column (integer or empty) if present."""
    if len(row) >= 3 and (row[0].strip() == "" or row[0].strip().isdigit()):
        return row[1:]
    return row


def _parse_rename_map(path: Optional[str]) -> Dict[str, str]:
    """Read a two-column (old_name, new_name) mapping; header optional."""
    out: Dict[str, str] = {}
    if not path:
        return out
    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            row = _strip_index(row)
            if len(row) < 2:
                continue
            key, val = row[0].strip(), row[1].strip()
            if not key or key.lower() in ("old_name", "gene_id", "loc_gene", "gene"):
                continue
            # accept both 'gene-NAME' and 'NAME' keys
            for k in (key, key.removeprefix("gene-")):
                out.setdefault(k, val)
    return out


def _parse_add_list(path: Optional[str]) -> Dict[str, str]:
    """Read a two-column (de-novo gene_id, new_name) list; header optional."""
    out: Dict[str, str] = {}
    if not path:
        return out
    with open(path, "r", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            row = _strip_index(row)
            if len(row) < 2:
                continue
            gid, name = row[0].strip(), row[1].strip()
            if not gid or gid.lower() in ("gene_id", "gene", "id"):
                continue
            out[gid] = name
    return out


def _parse_replace_list(path: Optional[str]) -> List[Tuple[str, str]]:
    """Read a two-column (liftoff_gene, new_name) replace list; header optional."""
    out: List[Tuple[str, str]] = []
    if not path:
        return out
    with open(path, "r", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            row = _strip_index(row)
            if len(row) < 2:
                continue
            lo, name = row[0].strip(), row[1].strip()
            if not lo or lo.lower() in ("loc_gene", "liftoff_gene", "gene", "old_name"):
                continue
            out.append((lo, name))
    return out


def _write_decisions(path: str, result) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["category", "liftoff_gene", "denovo_gene", "detail"])
        for d in result.decisions:
            writer.writerow([d.category, d.liftoff_gene, d.denovo_gene, d.detail])


def _write_categories(outdir: str, prefix: str, result) -> None:
    by: Dict[str, List] = {
        "replaced_loc": [], "replaced_concordant": [], "discordant_kept": [],
        "replace_skipped_discordant_kept": [], "added": [], "renamed": [],
    }
    for d in result.decisions:
        if d.category in by:
            by[d.category].append(d)
    name_map = {
        "replaced_loc": "replaced_by_de_novo_LOC_rule",
        "replaced_concordant": "replaced_by_de_novo_concordant",
        "discordant_kept": "discordant_kept",
        "replace_skipped_discordant_kept": "replace_skipped_discordant_kept",
        "added": "added_novel",
        "renamed": "renamed",
    }
    for cat, rows in by.items():
        if not rows:
            continue
        out = os.path.join(outdir, f"{prefix}_{name_map[cat]}.tsv")
        with open(out, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["category", "liftoff_gene", "denovo_gene", "detail"])
            for d in rows:
                writer.writerow([d.category, d.liftoff_gene, d.denovo_gene, d.detail])


def _write_report(path: str, result) -> None:
    counts = result.counts
    discordant = counts.get("discordant_kept", 0)
    loc = counts.get("replaced_loc", 0)
    conc = counts.get("replaced_concordant", 0)
    renamed = counts.get("renamed", 0)
    skip_disc = counts.get("replace_skipped_discordant_kept", 0)
    skip_nodn = counts.get("replace_skipped_no_denovo", 0)
    added = counts.get("added", 0)
    no_name = counts.get("denovo_skipped_no_name", 0)
    no_cds_len = counts.get("denovo_skipped_cds_length", 0)
    no_cds_ov = counts.get("denovo_skipped_cds_overlap", 0)
    dup_name = counts.get("denovo_skipped_dup_base_name", 0)
    no_cds = counts.get("denovo_skipped_no_cds", 0)
    no_expr = counts.get("denovo_skipped_expression", 0)
    lines = [
        "Merge report",
        "============",
        "",
        "Liftoff genes:",
        f"  replaced (LOC rule)           : {loc}",
        f"  replaced (concordant)         : {conc}",
        f"  renamed (name-change rule)    : {renamed}",
        f"  kept (discordant name, fix)   : {skip_disc}",
        f"  kept (no de-novo found)       : {skip_nodn}",
        "",
        "De-novo genes:",
        f"  added as novel                : {added}",
        f"  placed via replace            : {loc + conc}",
        f"  skipped (no name)             : {no_name}",
        f"  skipped (CDS length)          : {no_cds_len}",
        f"  skipped (CDS overlap)         : {no_cds_ov}",
        f"  skipped (dup base name)       : {dup_name}",
        f"  skipped (no CDS)              : {no_cds}",
        f"  skipped (expression)          : {no_expr}",
        "",
        f"Final genes written             : {len(result.final_genes)}",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--liftOff", required=True, help="Liftoff-transferred annotation (GFF3).")
    p.add_argument("--denovo", required=True,
                   help="De-novo annotation (GFF3) with functional names in the Name= attribute.")
    p.add_argument("--outdir", default="output", help="Directory for output files.")
    p.add_argument("--prefix", default="merged", help="Prefix for output filenames.")
    p.add_argument("--rename-map", default=None,
                   help="TSV (old_name<TAB>new_name) for the name-change rule.")
    p.add_argument("--add-list", default=None,
                   help="TSV (de-novo gene_id<TAB>new_name) of curated de-novo genes to add.")
    p.add_argument("--replace-list", default=None,
                   help="TSV (liftoff_gene<TAB>new_name) of curated Liftoff genes to replace.")
    p.add_argument("--list-mode", action="store_true",
                   help="Use the curated lists (rename-map/add-list/replace-list) "
                        "instead of automatic rule derivation.")
    p.add_argument("--require-name", action="store_true", default=True,
                   help="Only consider de-novo genes with a functional name (default).")
    p.add_argument("--no-require-name", dest="require_name", action="store_false",
                   help="Also consider de-novo genes without a functional name.")
    p.add_argument("--min-cds-bp", type=int, default=0,
                   help="Minimum de-novo CDS length (bp) to consider (default 0).")
    p.add_argument("--overlap-fraction", type=float, default=0.0,
                   help="Minimum gene-span overlap fraction to trigger the replace rule (default 0.0 = any).")
    p.add_argument("--cds-overlap-fraction", type=float, default=0.0,
                   help="Maximum CDS overlap fraction a novel gene may have with retained Liftoff genes (default 0.0 = none).")
    p.add_argument("--no-skip-dup-base-name", dest="skip_duplicate_base_name",
                   action="store_false", default=True,
                   help="Allow adding a de-novo gene whose base name already exists in the Liftoff set.")
    p.add_argument("--bigwig", default=None,
                   help="Optional BigWig expression track for the add rule (requires pyBigWig).")
    p.add_argument("--min-expression", type=float, default=0.0,
                   help="Minimum mean BigWig signal to keep a de-novo gene in the add rule.")
    p.add_argument("--bigwig-chrom-style", default="ucsc",
                   choices=["ucsc", "ensembl"],
                   help="BigWig chromosome naming convention (default ucsc).")
    p.add_argument("--version", action="version", version="annotation-merging 1.0.0")
    args = p.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[1/3] Parsing Liftoff annotation: {args.liftOff}", flush=True)
    liftoff = parse_gff(args.liftOff)
    print(f"      {len(liftoff)} genes", flush=True)

    print(f"[2/3] Parsing de-novo annotation: {args.denovo}", flush=True)
    denovo = parse_gff(args.denovo)
    print(f"      {len(denovo)} genes", flush=True)

    cfg = MergeConfig(
        require_name=args.require_name,
        min_cds_bp=args.min_cds_bp,
        overlap_fraction=args.overlap_fraction,
        cds_overlap_fraction=args.cds_overlap_fraction,
        skip_duplicate_base_name=args.skip_duplicate_base_name,
        rename_map=_parse_rename_map(args.rename_map),
        bigwig=args.bigwig,
        min_expression=args.min_expression,
        bigwig_chrom_style=args.bigwig_chrom_style,
    )

    use_lists = args.list_mode or args.add_list or args.replace_list

    print("[3/3] Merging (name-change -> replace -> add)...", flush=True)
    if use_lists:
        rename_map = _parse_rename_map(args.rename_map)
        add_list = _parse_add_list(args.add_list)
        replace_list = _parse_replace_list(args.replace_list)
        print(f"      list mode: rename={len(rename_map)} add={len(add_list)} "
              f"replace={len(replace_list)}", flush=True)
        result = merge_with_lists(
            liftoff, denovo, cfg,
            rename_map=rename_map, add_list=add_list, replace_list=replace_list)
    else:
        result = merge_annotations(liftoff, denovo, cfg)

    out_gff = os.path.join(args.outdir, f"{args.prefix}_merged.gff3")
    out_dec = os.path.join(args.outdir, f"{args.prefix}_decisions.tsv")
    out_rep = os.path.join(args.outdir, f"{args.prefix}_report.txt")

    write_gff(result.final_genes, out_gff,
              header=["source=annotation-merging",
                      "input_liftoff=" + os.path.basename(args.liftOff),
                      "input_denovo=" + os.path.basename(args.denovo)])
    _write_decisions(out_dec, result)
    _write_categories(args.outdir, args.prefix, result)
    _write_report(out_rep, result)

    print(f"      final genes: {len(result.final_genes)}", flush=True)
    print(f"      report:      {out_rep}", flush=True)
    print(f"      GFF3:        {out_gff}", flush=True)
    print(f"      decisions:   {out_dec}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
