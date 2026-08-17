#!/usr/bin/env python3
"""Re-derive the peaks2utr "changed-gene" subset from the merged GFF.

A gene is "changed" when its CDS span differs from the previous peaks2utr input
(the old annotation).  Only those genes need fresh 3' UTRs from peaks2utr; the
rest reuse the previous peaks2utr output's UTRs (see ``assemble_final.py``).

Run after the merge (which now normalizes Parent linkage) so the subset carries
correct ``gene -> mRNA -> CDS`` relationships, which peaks2utr's gffutils db
needs in order to find transcripts.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from src.gffio import parse_gff, write_gff


def cds_span(g):
    cds = [(r[3], r[4]) for r in g.child_rows if r[2] == "CDS"]
    if cds:
        return (g.seqid, min(s for s, e in cds), max(e for s, e in cds))
    exons = [(r[3], r[4]) for r in g.child_rows if r[2] == "exon"]
    if exons:
        return (g.seqid, min(s for s, e in exons), max(e for s, e in exons))
    return (g.seqid, g.start, g.end)


def main(old_input, merged, out):
    print("parsing old input...", flush=True)
    old_in = parse_gff(old_input)
    old_cds = {cds_span(g) for g in old_in}
    del old_in

    print("parsing merged GFF...", flush=True)
    new = parse_gff(merged)

    changed = [g for g in new if cds_span(g) not in old_cds]
    unchanged = [g for g in new if cds_span(g) in old_cds]
    print(f"changed: {len(changed)}  unchanged: {len(unchanged)}  total: {len(new)}", flush=True)

    write_gff(changed, out, header=["source=annotation-merging; peaks2utr subset of changed-CDS-span genes"])
    print(f"wrote {out}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--old-input", required=True,
                   help="Previous peaks2utr input GFF (the old annotation).")
    p.add_argument("--merged", required=True,
                   help="Merged annotation GFF.")
    p.add_argument("--out", required=True,
                   help="Output GFF of changed-CDS-span genes (for peaks2utr).")
    a = p.parse_args()
    main(a.old_input, a.merged, a.out)
