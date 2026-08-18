#!/usr/bin/env python
"""Concatenate trimmed per-gene alignments into a PHYLIP supermatrix (step 4).

Concatenates each species' trimmed alignment across all shared genes in the
order recorded by ``--gene-list``, producing a relaxed sequential PHYLIP file
suitable as RAxML input.

Usage (values come from config.sh):
    python 04_concat_supermatrix.py \
        --species-list "$EVO_SPECIES_LIST" \
        --gene-list    "$EVO_OUT_DIR/shared_singlecopy_genes.txt" \
        --trim-dir     "$EVO_OUT_DIR/03_alignments" \
        --out          "$EVO_OUT_DIR/04_supermatrix/supermatrix.phy"
"""

import argparse
import os
import sys


def parse_species_list(path):
    labels = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.lstrip().startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 2:
                labels.append(cols[1].strip())
    return labels


def read_fasta(path):
    """Return {header: sequence} for a FASTA file (headers keyed by first token)."""
    seqs = {}
    cur = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species-list", required=True, help="species_list.tsv")
    ap.add_argument("--gene-list", required=True, help="ordered shared-gene id list")
    ap.add_argument("--trim-dir", required=True, help="dir of <gene>.trim.faa alignments")
    ap.add_argument("--out", required=True, help="output relaxed-PHYLIP supermatrix")
    args = ap.parse_args()

    labels = parse_species_list(args.species_list)
    genes = [line.strip() for line in open(args.gene_list) if line.strip()]

    supermatrix = {label: [] for label in labels}
    total_len = 0
    for gene in genes:
        path = os.path.join(args.trim_dir, gene + ".trim.faa")
        if not os.path.isfile(path):
            sys.exit("ERROR: missing %s" % path)
        seqs = read_fasta(path)
        missing = [lbl for lbl in labels if lbl not in seqs]
        if missing:
            sys.exit("ERROR: %s missing from %s" % (", ".join(missing), path))
        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            sys.exit("ERROR: unequal sequence lengths in %s: %s" % (path, sorted(lengths)))
        total_len += lengths.pop()
        for lbl in labels:
            supermatrix[lbl].append(seqs[lbl])

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w") as out:
        out.write("%d %d\n" % (len(labels), total_len))
        for lbl in labels:
            out.write("%s %s\n" % (lbl, "".join(supermatrix[lbl])))

    print("Supermatrix: %d taxa x %d amino acids -> %s" % (len(labels), total_len, args.out))


if __name__ == "__main__":
    main()
