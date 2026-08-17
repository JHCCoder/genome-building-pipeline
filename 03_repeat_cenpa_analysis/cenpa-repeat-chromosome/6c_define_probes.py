#!/usr/bin/env python
"""
6c_define_probes.py -- define family-specific 31-mer probe sets.

A probe is a 31-mer that is diagnostic of ONE repeat family:
  * present in the family's array sequence (count >= C_FAM)
  * absent from the OTHER repeat families' array sequences (count == 0)
  * enriched in the family relative to the whole genome (genome_count / family
    count not excessive) -- this guards against k-mers that are common genome
    motifs rather than family-specific sequence.
This makes the probe "sequence-specific to the family".

Outputs (per family): a file of 31-mers (one per line), capped at MAX_PROBES
ranked by family count, plus a small fasta for KMC probe-DB construction.

NOTE: KMC counters are capped at 255 by default (-cs); both family and genome
counts share the same cap, so the ratio filter uses raw capped values which is
conservative (it can only miss, never wrongly include, high-copy k-mers).

Usage: python 6c_define_probes.py <probe_dir>
"""
import sys, os

C_FAM = 5          # min copies of the k-mer within the family arrays
MAX_PROBES = 5000  # cap probe set size
MAX_GENOME_RATIO = 20  # genome_count/family_count <= this (not a common motif)

def load_counts(path):
    d = {}
    if not os.path.exists(path):
        return d
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == 2:
                d[p[0]] = int(p[1])
    return d

def probes_for(family_counts, other_counts_list, genome_counts, out_prefix):
    others = [set(c.keys()) for c in other_counts_list]
    union_other = set().union(*others) if others else set()
    cand = []
    for k, v in family_counts.items():
        if v < C_FAM or k in union_other:
            continue
        g = genome_counts.get(k)
        if g is not None and g / max(v, 1) > MAX_GENOME_RATIO:
            continue   # k-mer is a common genome motif, not family-specific
        cand.append((k, v))
    cand.sort(key=lambda x: -x[1])
    cand = cand[:MAX_PROBES]
    with open(out_prefix + ".31mers.txt", "w") as f:
        for k, v in cand:
            f.write(k + "\n")
    with open(out_prefix + ".fa", "w") as f:
        for i, (k, v) in enumerate(cand):
            f.write(f">probe{i}\n{k}\n")
    return len(cand)

def main():
    probe_dir = sys.argv[1] if len(sys.argv) > 1 else \
        "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome/data/probes"

    b6 = load_counts(os.path.join(probe_dir, "bin6_counts.tsv"))   # 349
    b4 = load_counts(os.path.join(probe_dir, "bin4_counts.tsv"))   # 195
    b8 = load_counts(os.path.join(probe_dir, "bin8_counts.tsv"))   # 389
    g6 = load_counts(os.path.join(probe_dir, "bin6_genomecounts.tsv"))
    g4 = load_counts(os.path.join(probe_dir, "bin4_genomecounts.tsv"))
    g8 = load_counts(os.path.join(probe_dir, "bin8_genomecounts.tsv"))

    n349 = probes_for(b6, [b4, b8], g6, os.path.join(probe_dir, "probes_349"))
    n195 = probes_for(b4, [b6, b8], g4, os.path.join(probe_dir, "probes_195"))
    n389 = probes_for(b8, [b4, b6], g8, os.path.join(probe_dir, "probes_389"))
    print(f"probes_349: {n349}")
    print(f"probes_195: {n195}")
    print(f"probes_389: {n389}")

    # L1 subfamily probes from consensi (fam189 vs fam18).
    # A consensus is a single sequence, so each 31-mer occurs 1-4x; use
    # count >= 2 (i.e. present in >=2 tandem copies) and absent from the
    # other L1 subfamily consensus.
    C_CONS = 2
    fam189 = load_counts(os.path.join(probe_dir, "fam189_counts.tsv"))
    fam18 = load_counts(os.path.join(probe_dir, "fam18_counts.tsv"))
    with open(os.path.join(probe_dir, "probes_fam189.31mers.txt"), "w") as f:
        for k in [k for k in fam189 if fam189[k] >= C_CONS and k not in fam18][:MAX_PROBES]:
            f.write(k + "\n")
    with open(os.path.join(probe_dir, "probes_fam18.31mers.txt"), "w") as f:
        for k in [k for k in fam18 if fam18[k] >= C_CONS and k not in fam189][:MAX_PROBES]:
            f.write(k + "\n")
    print("fam189/fam18 subfamily probes written (C_CONS=%d)" % C_CONS)

if __name__ == "__main__":
    main()
