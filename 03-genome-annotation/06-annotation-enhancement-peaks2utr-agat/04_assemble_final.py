#!/usr/bin/env python3
"""Assemble the final annotation: merged GFF + peaks2utr 5'/3' UTRs.

For each gene in the new annotation, the peaks2utr UTR features are taken from
the run that matches its model:
  * unchanged genes (CDS span identical to the previous peaks2utr input) use the
    previous peaks2utr output;
  * changed genes (the subset) use the new peaks2utr subset output.

UTR features are re-keyed to the new annotation's transcript IDs (Liftoff
`rna-XM_*` IDs already match; de-novo `gene-<name>.t<n>` IDs are matched to the
new Braker transcript via their CDS span).
"""
import os
import re
import subprocess
import sys
from collections import defaultdict

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
OLD_OUT = _cfg("PEAKS2UTR_PREV_OUTPUT")
NEW = _cfg("MERGED_GFF")
NEW_OUT = _cfg("PEAKS2UTR_SUBSET_OUT")
FINAL = _cfg("FINAL_UTR_GFF")

# The GFF3 helpers live in the sibling annotation-merging tool (step 05).
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "05-annotation-merging"))
from src.gffio import parse_gff, write_gff, set_attributes

UTR_TYPES = {"five_prime_UTR", "three_prime_UTR"}
TRANSCRIPT_TYPES = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA", "snRNA", "snoRNA"}


def cds_span(g):
    cds = [(r[3], r[4]) for r in g.child_rows if r[2] == "CDS"]
    if cds:
        return (g.seqid, min(s for s, e in cds), max(e for s, e in cds))
    exons = [(r[3], r[4]) for r in g.child_rows if r[2] == "exon"]
    if exons:
        return (g.seqid, min(s for s, e in exons), max(e for s, e in exons))
    return (g.seqid, g.start, g.end)


def transcript_cds_map(g):
    """{mRNA_id: cds_span} for a gene, from CDS rows grouped by Parent."""
    out = {}
    cds_by_tx = defaultdict(list)
    for r in g.child_rows:
        if r[2] == "CDS":
            m = re.search(r"Parent=([^;]+)", r[8])
            if m:
                cds_by_tx[m.group(1)].append((r[3], r[4]))
    for tx, cds in cds_by_tx.items():
        out[tx] = (min(s for s, e in cds), max(e for s, e in cds))
    return out


def build_source_index(source_genes):
    """Map source (peaks2utr-output) genes by CDS span, with per-gene
    {transcript_id: cds_span} and UTR features."""
    by_cds = {}
    for g in source_genes:
        by_cds[cds_span(g)] = g
    return by_cds


def transcript_index(g):
    """{transcript_id: (start, end, strand)} for every transcript-type feature
    of a gene.  Unlike :func:`transcript_cds_map`, this does not depend on CDS
    ``Parent`` attributes (de-novo genes from ``braker_peak2utr.gff3`` have
    orphaned CDS/exon rows with no Parent, and ncRNA has no CDS at all)."""
    out = {}
    for r in g.child_rows:
        if r[2] in TRANSCRIPT_TYPES:
            m = re.search(r"ID=([^;]+)", r[8])
            if m:
                out[m.group(1)] = (r[3], r[4], r[6])
    return out


def resolve_transcript(new_gene, old_tx, old_tx_map, new_tx_map, utr_row):
    """Resolve which transcript of *new_gene* a UTR (old ``Parent=old_tx``)
    belongs to, tried in order of decreasing confidence:

      1. direct ID match — Liftoff/RefSeq transcript IDs are stable;
      2. CDS span match — for coding transcripts whose CDS carries a Parent;
      3. single-transcript fallback — the common de-novo / ncRNA case;
      4. nearest 3' end — positional last resort for multi-transcript genes.
    """
    tx = transcript_index(new_gene)
    if old_tx in tx:
        return old_tx
    old_cds = old_tx_map.get(old_tx)
    if old_cds:
        for t, span in new_tx_map.items():
            if span == old_cds:
                return t
        for t, span in new_tx_map.items():
            if span[0] <= old_cds[1] and span[1] >= old_cds[0]:
                return t
    if len(tx) == 1:
        return next(iter(tx))
    # positional: closest transcript 3' end to the UTR midpoint
    u_mid = (utr_row[3] + utr_row[4]) // 2
    best, best_d = None, None
    for t, (s, e, strand) in tx.items():
        three_end = e if strand == "+" else s
        d = abs(three_end - u_mid)
        if best_d is None or d < best_d:
            best, best_d = t, d
    return best


def main():
    print("parsing old input...", flush=True)
    old_in = parse_gff(OLD_IN)
    old_in_cds = {cds_span(g) for g in old_in}
    del old_in

    print("parsing old peaks2utr output...", flush=True)
    old_out = parse_gff(OLD_OUT)
    old_out_by_cds = build_source_index(old_out)

    print("parsing new GFF...", flush=True)
    new = parse_gff(NEW)

    # optional: new subset peaks2utr output (when it exists)
    new_out_by_cds = None
    new_out_by_name = None
    if os.path.isfile(NEW_OUT):
        print("parsing new peaks2utr subset output...", flush=True)
        new_out = parse_gff(NEW_OUT)
        new_out_by_cds = build_source_index(new_out)
        # peaks2utr/gffutils keeps only the outermost transcript per gene, so the
        # subset output's gene-level CDS span is narrower than the merged gene's
        # (multi-transcript) span.  Match changed genes by gene ID instead, which
        # the subset output preserves unchanged.
        new_out_by_name = {g.gene_id.removeprefix("gene-").lower(): g for g in new_out}
    else:
        print(f"note: {NEW_OUT} not ready; changed genes will keep no UTRs", flush=True)

    final = []
    n_unchanged = n_changed = n_utrs = n_skipped = 0
    for g in new:
        cds = cds_span(g)
        changed = cds not in old_in_cds
        if changed and new_out_by_name is None:
            # subset peaks2utr output not ready; changed genes keep no UTRs
            final.append(g)
            n_changed += 1
            continue
        if changed:
            source_gene = new_out_by_name.get(g.gene_id.removeprefix("gene-").lower())
        else:
            source_gene = old_out_by_cds.get(cds)
        if source_gene is None:
            n_skipped += 1
            final.append(g)
            continue
        new_tx_map = transcript_cds_map(g)
        old_tx_map = transcript_cds_map(source_gene)
        for r in source_gene.child_rows:
            if r[2] not in UTR_TYPES:
                continue
            m = re.search(r"Parent=([^;]+)", r[8])
            if m:
                old_tx = m.group(1)
                new_tx = resolve_transcript(g, old_tx, old_tx_map, new_tx_map, r)
                if new_tx:
                    attrs = set_attributes(r[8], Parent=new_tx)
                    g.child_rows.append((r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], attrs))
                    n_utrs += 1
                    continue
            # couldn't re-key: append with original Parent (may be dangling) as last resort
            g.child_rows.append(r)
            n_utrs += 1
        final.append(g)
        if changed:
            n_changed += 1
        else:
            n_unchanged += 1

    write_gff(final, FINAL, header=["source=annotation-merging + peaks2utr; final degu annotation"])
    print(f"\nFinal: {len(final)} genes; unchanged w/ UTRs: {n_unchanged}, changed: {n_changed}, "
          f"UTR features added: {n_utrs}, source-miss: {n_skipped}")
    print(f"wrote {FINAL}")


if __name__ == "__main__":
    main()
