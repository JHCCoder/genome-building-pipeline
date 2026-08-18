#!/usr/bin/env python3
"""Convert ``mash triangle`` output (lower-triangular) to a square PHYLIP matrix.

``mash triangle`` writes the pairwise Jaccard distances as a lower-triangular
matrix: a first line with the number of taxa, then one line per taxon holding
its name followed by the distances to all *previous* taxa.  FastME expects a
square PHYLIP matrix (full symmetric matrix with a leading taxon count), which
this script produces.

Usage:
    python3 convert_to_square.py all_distances_cleaned.tab all_distances_cleaned.phy
"""

import sys

import numpy as np


def convert_lower_tri_to_square(input_file, output_file):
    """Convert lower triangular distance matrix to square PHYLIP format."""
    with open(input_file) as f:
        lines = [line.strip() for line in f if line.strip()]

    # First line is the number of taxa.
    n_taxa = int(lines[0])
    names = []
    distances = []

    # Parse the lower triangular matrix.
    for line in lines[1:n_taxa + 1]:
        parts = line.split()
        names.append(parts[0])
        distances.append(list(map(float, parts[1:])))

    # Build the square matrix, mirroring to the upper triangle.
    matrix = np.zeros((n_taxa, n_taxa))
    for i in range(n_taxa):
        for j in range(i):
            matrix[i][j] = distances[i][j]
            matrix[j][i] = distances[i][j]

    with open(output_file, "w") as f:
        f.write(f"{n_taxa}\n")
        for i in range(n_taxa):
            f.write(f"{names[i]}")
            for j in range(n_taxa):
                f.write(f" {matrix[i][j]:.6f}")
            f.write("\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_to_square.py input_file output_file")
        sys.exit(1)
    convert_lower_tri_to_square(sys.argv[1], sys.argv[2])
    print(f"Successfully converted {sys.argv[1]} to square matrix format in {sys.argv[2]}")
