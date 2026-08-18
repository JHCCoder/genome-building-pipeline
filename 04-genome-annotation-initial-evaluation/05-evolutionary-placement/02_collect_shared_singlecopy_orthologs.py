#!/usr/bin/env python
"""Collect the single-copy orthologs shared by every species (step 2).

Reads each species' BUSCO ``full_table.tsv``, keeps the genes classified
"Complete" (BUSCO's single-copy complete class) in *every* genome, and writes
one multi-species FASTA per retained gene. Each sequence is relabelled to its
species ``tree_label`` so the downstream MAFFT -> trimAl -> concatenation steps
can key on a stable name.

Usage (values come from config.sh):
    python 02_collect_shared_singlecopy_orthologs.py \
        --species-list "$EVO_SPECIES_LIST" \
        --busco-dir    "$EVO_OUT_DIR/01_busco" \
        --lineage      "$EVO_BUSCO_LINEAGE" \
        --out-dir      "$EVO_OUT_DIR/02_single_copy_orthologs" \
        --gene-list    "$EVO_OUT_DIR/shared_singlecopy_genes.txt"
"""

import argparse
import os
import sys


def parse_species_list(path):
    """Return the ordered list of tree labels from species_list.tsv."""
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


def parse_complete_busco_ids(full_table):
    """Return the set of Busco ids with Status == "Complete" in a full_table.tsv."""
    complete = set()
    with open(full_table) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            # full_table.tsv columns: Busco id, Status, Sequence, Gene Start, ...
            if len(cols) > 1 and cols[1] == "Complete":
                complete.add(cols[0])
    return complete


def read_first_sequence(path):
    """Return the first FASTA record's sequence as one unwrapped string."""
    seq = []
    seen_header = False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header:
                    break
                seen_header = True
                continue
            seq.append(line)
    return "".join(seq)


def wrap(seq, width=60):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--species-list", required=True,
                    help="species_list.tsv (genome<tab>label)")
    ap.add_argument("--busco-dir", required=True,
                    help="dir containing one <label>/run_<lineage>/ per species")
    ap.add_argument("--lineage", required=True, help="BUSCO lineage dataset name")
    ap.add_argument("--out-dir", required=True, help="dir for per-gene multi-species FASTAs")
    ap.add_argument("--gene-list", required=True, help="output: ordered shared-gene id list")
    args = ap.parse_args()

    labels = parse_species_list(args.species_list)
    if not labels:
        sys.exit("ERROR: no species found in %s" % args.species_list)

    complete_per_species = {}
    seq_dir_per_species = {}
    for label in labels:
        run_dir = os.path.join(args.busco_dir, label, "run_%s" % args.lineage)
        full_table = os.path.join(run_dir, "full_table.tsv")
        if not os.path.isfile(full_table):
            sys.exit("ERROR: missing %s (did step 1 finish for %s?)" % (full_table, label))
        complete_per_species[label] = parse_complete_busco_ids(full_table)
        seq_dir_per_species[label] = os.path.join(
            run_dir, "busco_sequences", "single_copy_busco_sequences"
        )

    shared = set(complete_per_species[labels[0]])
    for label in labels[1:]:
        shared &= complete_per_species[label]
    shared = sorted(shared)

    os.makedirs(args.out_dir, exist_ok=True)
    written = 0
    for gene in shared:
        out_fa = os.path.join(args.out_dir, gene + ".faa")
        with open(out_fa, "w") as out:
            for label in labels:
                seq_path = os.path.join(seq_dir_per_species[label], gene + ".faa")
                if not os.path.isfile(seq_path):
                    sys.exit("ERROR: missing %s" % seq_path)
                out.write(">%s\n%s\n" % (label, wrap(read_first_sequence(seq_path))))
        written += 1

    with open(args.gene_list, "w") as g:
        for gene in shared:
            g.write(gene + "\n")

    print("%d species, %d shared single-copy orthologs" % (len(labels), written))
    print("per-gene FASTAs -> %s" % args.out_dir)
    print("gene order      -> %s" % args.gene_list)


if __name__ == "__main__":
    main()
