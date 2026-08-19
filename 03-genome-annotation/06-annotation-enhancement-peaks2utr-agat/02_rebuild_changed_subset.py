#!/usr/bin/env python3
"""Re-derive the peaks2utr "changed-gene" subset from the merged GFF.

A gene is "changed" when its CDS span differs from the previous peaks2utr input
(``PEAKS2UTR_PREV_INPUT``, i.e. the old annotation).  Only those genes need
fresh 3' UTRs from peaks2utr; the rest reuse the previous peaks2utr output's
UTRs (see ``04_assemble_final.py``).

Run after the merge (which normalizes Parent linkage) so the subset carries
correct ``gene -> mRNA -> CDS`` relationships, which peaks2utr's gffutils db
needs in order to find transcripts.
"""
import os
import subprocess
import sys

# --- Paths come from config.sh at the repo root (override via the environment).
def _cfg(name):
    val = os.environ.get(name)
    if val:
        return val
    root = os.path.abspath(os.path.dirname(__file__))
    while not os.path.isfile(os.path.join(root, "config.sh")) and root != os.path.dirname(root):
        root = os.path.dirname(root)
    # Source config.sh in a subshell so nested $VAR references are expanded.
    val = subprocess.check_output(
        ["bash", "-c", f'source "$1" 2>/dev/null; printf "%s" "${{{name}}}"', "_",
         os.path.join(root, "config.sh")],
        text=True,
    )
    if not val:
        sys.exit(f"ERROR: {name} is not set in the environment or in config.sh")
    return val


OLD_IN = _cfg("PEAKS2UTR_PREV_INPUT")
NEW = _cfg("MERGED_GFF")
OUT = _cfg("CHANGED_SUBSET_GFF")

# The GFF3 helpers live in the sibling annotation-merging tool (step 05).
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "05-annotation-merging"))
from src.gffio import parse_gff, write_gff


def cds_span(g):
    cds = [(r[3], r[4]) for r in g.child_rows if r[2] == "CDS"]
    if cds:
        return (g.seqid, min(s for s, e in cds), max(e for s, e in cds))
    exons = [(r[3], r[4]) for r in g.child_rows if r[2] == "exon"]
    if exons:
        return (g.seqid, min(s for s, e in exons), max(e for s, e in exons))
    return (g.seqid, g.start, g.end)


def main():
    print("parsing old input...", flush=True)
    old_in = parse_gff(OLD_IN)
    old_cds = {cds_span(g) for g in old_in}
    del old_in

    print("parsing merged GFF...", flush=True)
    new = parse_gff(NEW)

    changed = [g for g in new if cds_span(g) not in old_cds]
    unchanged = [g for g in new if cds_span(g) in old_cds]
    print(f"changed: {len(changed)}  unchanged: {len(unchanged)}  total: {len(new)}", flush=True)

    write_gff(changed, OUT, header=["source=annotation-merging; peaks2utr subset of changed-CDS-span genes"])
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
