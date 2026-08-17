#!/usr/bin/env python
"""Compute basic assembly metrics for a FASTA file.

Extracted from code/jupyter-notebook-script/compare_fasta_files.ipynb.

Metrics: number of bases, number of sequences, N50, L50, largest scaffold,
GC content (total and ungapped), and the count / % of sequence length above
the 25 Mbp / 10 Mbp / 1 Mbp / 15 kbp thresholds.

Requires Biopython (Bio.SeqIO).

Usage:
    python 04_fasta_metrics.py <fasta_file>
"""

import sys
from Bio import SeqIO


def count_bases(fasta_file):
    length = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            else:
                length += len(line.strip())
    return length


def count_sequences(fasta_file):
    count = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                count += 1
    return count


def calculate_N50(fasta_file):
    sequence_lengths = []
    with open(fasta_file, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequence_lengths.append(len(sequence))
                    sequence = ''
            else:
                sequence += line.strip()
        if sequence:  # Add the last sequence
            sequence_lengths.append(len(sequence))

    sequence_lengths.sort(reverse=True)
    total_length = sum(sequence_lengths)

    cumulative_length = 0
    l50 = 0
    for length in sequence_lengths:
        cumulative_length += length
        l50 += 1
        if cumulative_length >= total_length / 2:
            return length, sequence_lengths, l50


def calculate_total_gc_content(fasta_file):
    """Average GC content across all sequences (gaps included)."""
    total_gc = 0
    total_length = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        gc_count = sequence.count('G') + sequence.count('C')
        total_gc += gc_count
        total_length += len(sequence)
    if total_length == 0:
        return 0.0
    return (total_gc / total_length) * 100


def calculate_ungapped_gc_content(fasta_file):
    """Average GC content ignoring N / - gap characters."""
    total_gc = 0
    total_ungapped_length = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        ungapped_seq = seq.replace('N', '').replace('-', '')
        gc_count = ungapped_seq.count('G') + ungapped_seq.count('C')
        total_gc += gc_count
        total_ungapped_length += len(ungapped_seq)
    if total_ungapped_length == 0:
        return 0.0
    return (total_gc / total_ungapped_length) * 100


def print_metrics(fasta_file):
    bp_num = count_bases(fasta_file)
    print("Number of base pairs in the FASTA file is:", bp_num)

    seq_num = count_sequences(fasta_file)
    print("Number of sequences in the FASTA file is:", seq_num)

    n50 = calculate_N50(fasta_file)
    print("N50 of the FASTA file is:", n50[0])
    print("L50 of the FASTA file is:", n50[2])
    print("Largest scaffold is:", n50[1][0])

    gc_results = calculate_ungapped_gc_content(fasta_file)
    print("GC content of the FASTA file is:", gc_results)

    l1 = [i for i in n50[1] if i > 25000000]
    perc_25mbp = sum(l1) / bp_num
    print("# of sequence > 25Mbp:", len(l1))
    print("% of sum length of sequence >25Mbp:", perc_25mbp)

    l2 = [i for i in n50[1] if i > 10000000]
    perc_10mbp = sum(l2) / bp_num
    print("# of sequence > 10Mbp:", len(l2))
    print("% of sum length of sequence >10Mbp:", perc_10mbp)

    l3 = [i for i in n50[1] if i > 1000000]
    perc_1mbp = sum(l3) / bp_num
    print("# of sequence > 1Mbp:", len(l3))
    print("% of sum length of sequence > 1Mbp:", perc_1mbp)

    l4 = [i for i in n50[1] if i > 15000]
    perc_15kbp = sum(l4) / bp_num
    print("# of sequence > 15Kbp:", len(l4))
    print("% of sum length of sequence > 15Kbp:", perc_15kbp)
    print("##############")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: python 04_fasta_metrics.py <fasta_file>")
    print_metrics(sys.argv[1])
