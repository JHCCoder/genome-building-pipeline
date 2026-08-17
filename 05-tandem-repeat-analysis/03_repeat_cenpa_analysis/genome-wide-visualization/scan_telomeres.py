#!/usr/bin/env python3
"""
Scan the degu genome for telomeric tandem repeats.

Probes:
  - TTAGGG  (canonical vertebrate telomere, plus-strand)
  - CCCTAA  (reverse-complement orientation of TTAGGG arrays)
  - TTCAGGG / CCCTGAA (murid variant, tested for completeness)

For each target chromosome (chr1-28, chrX, chrY) we stream the assembly
FASTA, regex-scan the full sequence, and merge same-strand runs into
arrays. Outputs a BED-like TSV with per-array coordinates, strand,
length, and repeat count.

Run:  python3 scan_telomeres.py <assembly.fasta> <assembly.fai> \
            --chrom chr1 ... --min-copies 6 --out telomere_arrays.tsv
"""
import argparse
import re
import sys

# Probes: motif string -> regex compiled per min_copies in scan_chrom
PROBES = ["TTAGGG", "CCCTAA", "TTCAGGG", "CCCTGAA"]


def parse_fai(fai_path):
    """Return dict chrom -> length."""
    out = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out[parts[0]] = int(parts[1])
    return out


def stream_fasta(fasta_path, chroms):
    """Yield (chrom, sequence) for each requested chromosome, in file order."""
    wanted = set(chroms)
    cur_chrom = None
    cur_seq = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_chrom is not None and cur_chrom in wanted:
                    yield cur_chrom, "".join(cur_seq).upper()
                name = line[1:].split()[0]
                cur_chrom = name
                cur_seq = []
            else:
                cur_seq.append(line.rstrip("\n"))
    if cur_chrom is not None and cur_chrom in wanted:
        yield cur_chrom, "".join(cur_seq).upper()


def scan_chrom(chrom, seq, min_copies):
    arrays = []
    for motif in PROBES:
        regex = re.compile(r"(%s){%d,}" % (motif, min_copies))
        for m in regex.finditer(seq):
            start = m.start()
            end = m.end()
            n_copies = (end - start) / len(motif)
            arrays.append((chrom, start, end, motif, round(n_copies, 2)))
    arrays.sort(key=lambda x: (x[1], x[2]))
    return arrays


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("fai")
    ap.add_argument("--chrom", nargs="+", required=True)
    ap.add_argument("--min-copies", type=int, default=6)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    chroms = args.chrom
    min_copies = args.min_copies

    all_arrays = []
    chrom_lens = parse_fai(args.fai)

    for chrom, seq in stream_fasta(args.fasta, chroms):
        arrays = scan_chrom(chrom, seq, min_copies)
        # Annotate with chromosome length
        clen = chrom_lens.get(chrom, None)
        for a in arrays:
            # position relative to chromosome end (0 = at start, 1 = at end)
            dist_from_start = a[1]
            dist_from_end = (clen - a[2]) if clen else None
            all_arrays.append((a[0], a[1], a[2], a[3], a[4], clen, dist_from_start, dist_from_end))
        sys.stderr.write("done %s: %d arrays, seq len %d\n" % (chrom, len(arrays), len(seq)))

    with open(args.out, "w") as f:
        f.write("chrom\tstart\tend\tstrand_motif\tn_copies\tchrom_len\tdist_from_start\tdist_from_end\n")
        for r in all_arrays:
            f.write("\t".join(str(x) for x in r) + "\n")

    print("Total arrays (min_copies=%d): %d" % (min_copies, len(all_arrays)))


if __name__ == "__main__":
    main()
