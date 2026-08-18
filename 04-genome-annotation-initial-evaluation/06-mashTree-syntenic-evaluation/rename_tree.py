#!/usr/bin/env python3
"""Rename the leaf labels of a Newick tree using a two-column mapping file.

mashtree / mash / FastME write leaf labels using the raw input filenames
(e.g. ``GCF_000001635.27_GRCm39_genomic.fna``).  This script relabels them
to human-readable species names using a TSV with a header::

    old_name<TAB>new_name

Tree files produced by mashtree drop the ``.fna``/``.fasta`` extension, so each
mapping key is matched both with and without those extensions.  Replacements
are applied longest-key-first so that a key which is a prefix of another is not
partially rewritten.

Usage:
    python3 rename_tree.py input.nwk name_conversion.tsv output.nwk
"""

import re
import sys


def rename_tree_nodes(tree_file, mapping_file, output_file):
    # Read mapping file (skip the header line).
    name_map = {}
    with open(mapping_file) as f:
        next(f)
        for line in f:
            old, new = line.rstrip("\n").split("\t")
            name_map[old] = new
            # mashtree tree leaves have no .fna/.fasta extension.
            old_base = old.replace(".fna", "").replace(".fasta", "")
            name_map[old_base] = new

    with open(tree_file) as f:
        tree = f.read().strip()

    # Longest keys first so prefix keys don't clobber longer names.
    for old_name in sorted(name_map, key=len, reverse=True):
        new_name = name_map[old_name]
        # Match the name either immediately before a ':' (branch length) or ')'
        # (end of a leaf), or at word boundaries elsewhere.
        tree = re.sub(
            r"({})(?=[:)])|\b{}\b".format(re.escape(old_name), re.escape(old_name)),
            new_name,
            tree,
        )

    with open(output_file, "w") as f:
        f.write(tree)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rename_tree.py input_nwk name_conversion_file output_nwk")
        sys.exit(1)
    rename_tree_nodes(sys.argv[1], sys.argv[2], sys.argv[3])
    print(f"Renamed tree written to {sys.argv[3]}")
