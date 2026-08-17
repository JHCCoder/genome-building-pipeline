#!/usr/bin/env python3
"""Classify the 3,367 period-110-135 TRF hits against canonical minor/major
mouse satellite consensus.

KEY FIX (was wrong in earlier versions): GRCm39 genomic.fna is FASTA-wrapped
80 bp/line with a trailing newline (FAI lb=80, lw=81). The FASTA reader must
compute byte offsets as off + line*lw (NOT off + line*(lw+1)); the off-by-one-
per-line bug drifted large-coordinate extraction by ~1.2 Mb and corrupted every
earlier classification. Extraction is also case-insensitive (the fna is
soft-masked, so genuine satellite arrays are lowercase).

Discrimination (calibrated on correct extractions):
  minor:   SW vs canonical minor coverage >= 0.6 and > major coverage
  major:   SW vs canonical major (234bp GSAT_MM dimer) coverage >= 0.6 and > minor
  ambiguous: both >= 0.6, |cmi - cma| < 0.05 (minor/major share a GAAAC motif;
             the GSAT_MM dimer contains a minor-like half, so genuine minor can
             also align to the major canonical)
  other:   neither coverage >= 0.6 (TRF period false positive / unrelated repeat)

A score floor of 60 avoids calling highly divergent chance alignments satellite.
Score-only SW runs first; full-traceback (identity+coverage) SW runs only for
hits clearing the floor (parallel chunked via --start/--end/--out).
"""
import csv
import re
import sys

csv.field_size_limit(sys.maxsize)

CWD = "/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/code/command-line-script/genome-annotation/trf-tandem-repeat"
TSV = CWD + "/mouse_trf_minor120.tmp.tsv"
FA = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/GRCm39_genome/GCF_000001635.27_GRCm39_genomic.fna"
FAI = FA + ".fai"
OUT = CWD + "/classifier_3367_caseinsensitive.tsv"

MINOR = "GGAAAATGATAAAACAGAACTGTAGAACATATTAGATGAGTTCAGTTACAACTAAAAAACACATTGCTTGGAACCGCATTTTGTAGAACAGTGTATATCATGAGTTACAATGAGAAACAT"
MAJOR = "CCTGGAATATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTAGGACGTGAAATATGGCGAGGAAAACTGAAAAAGGTGGAAAATTTAGAAATGTCCACTGTAGGACGTGGAATATGGCAAGAAAACTGAAAATCATGGAAAATGAGAAACATCCACTTGACGACTTGAAAAATGACGAAATCACTAAAAAACGTGAAAAATGAGAAATGCACACTGAAGGA"

SCORE_FLOOR = 80  # applied to the WINNING canon's score (genuine satellite >= ~100)
COV_MIN = 0.6
COV_ZONE = 0.05


def load_fai(p):
    d = {}
    for line in open(p):
        f = line.split()
        d[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return d


def load_chr_map(fasta):
    m = {}
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            acc = line[1:].split()[0]
            mm = re.search(r"chromosome (\d+|X|Y),", line)
            if mm:
                m[f"chr{mm.group(1)}"] = acc
    return m


class FaiReader:
    """Reads [s:e) via the .fai. FAI columns: lb=bases/line, lw=bytes/line
    including the newline. Byte offset of line L is off + L*lw."""

    def __init__(self, path, fai, chr_map):
        self.fh = open(path)
        self.fai = fai
        self.chr_map = chr_map

    def get(self, name, s, e):
        acc = self.chr_map.get(name, name)
        L, off, lb, lw = self.fai[acc]
        line0 = s // lb
        base_in_line = s % lb
        bstart = off + line0 * lw
        nlines = (e - line0 * lb + lb - 1) // lb
        self.fh.seek(bstart)
        raw = self.fh.read(nlines * lw)
        seq = "".join(raw.split("\n"))
        return seq[base_in_line: base_in_line + (e - s)].upper()


def sw_score(seq, canon):
    """Score-only Smith-Waterman local (match=2, mismatch=-3, gap=-2)."""
    m, n = len(canon), len(seq)
    prev = [0] * (n + 1)
    best = 0
    for i in range(1, m + 1):
        cur = [0] * (n + 1)
        pj = prev
        cj = cur
        ci = canon[i - 1]
        for j in range(1, n + 1):
            d = pj[j - 1] + (2 if ci == seq[j - 1] else -3)
            u = pj[j] - 2
            l = cj[j - 1] - 2
            v = 0
            if d > v: v = d
            if u > v: v = u
            if l > v: v = l
            cj[j] = v
            if v > best:
                best = v
        prev = cur
    return best


def sw_full(seq, canon):
    """Full SW with traceback. Returns (score, identity, coverage)."""
    m, n = len(canon), len(seq)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    best = 0
    bi = bj = 0
    for i in range(1, m + 1):
        Hi = H[i]
        Hp = H[i - 1]
        ci = canon[i - 1]
        for j in range(1, n + 1):
            v = 0
            d = Hp[j - 1] + (2 if ci == seq[j - 1] else -3)
            u = Hp[j] - 2
            l = Hi[j - 1] - 2
            if d > v: v = d
            if u > v: v = u
            if l > v: v = l
            Hi[j] = v
            if v > best:
                best = v
                bi, bj = i, j
    i, j = bi, bj
    matches = 0
    pairs = 0
    while i > 0 and j > 0 and H[i][j] > 0:
        cur = H[i][j]
        if cur == H[i - 1][j - 1] + (2 if canon[i - 1] == seq[j - 1] else -3):
            pairs += 1
            if canon[i - 1] == seq[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif cur == H[i - 1][j] - 2:
            i -= 1
        elif cur == H[i][j - 1] - 2:
            j -= 1
        else:
            break
    cov = pairs / m if m else 0.0
    identity = matches / pairs if pairs else 0.0
    return best, identity, cov


def revcomp(s):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(c, "N") for c in reversed(s))


def classify(seq):
    rc = revcomp(seq)
    s1 = max(sw_score(seq, MINOR), sw_score(rc, MINOR))
    s2 = max(sw_score(seq, MAJOR), sw_score(rc, MAJOR))
    if max(s1, s2) < SCORE_FLOOR:
        return "other", s1, 0.0, 0.0, s2, 0.0, 0.0
    smi, imi, cmi = max((sw_full(seq, MINOR), sw_full(rc, MINOR)), key=lambda t: t[0])
    sma, ima, cma = max((sw_full(seq, MAJOR), sw_full(rc, MAJOR)), key=lambda t: t[0])
    if abs(cmi - cma) < COV_ZONE and max(cmi, cma) >= COV_MIN:
        # matches both canons substantially: genuine satellite unless weak
        return ("ambiguous" if max(smi, sma) >= SCORE_FLOOR else "other",
                smi, imi, cmi, sma, ima, cma)
    if cmi > cma:
        winner, wcov, wsc = "minor", cmi, smi
    else:
        winner, wcov, wsc = "major", cma, sma
    cat = winner if (wcov >= COV_MIN and wsc >= SCORE_FLOOR) else "other"
    return cat, smi, imi, cmi, sma, ima, cma


def main():
    a = sys.argv[1:]
    start = int(a[a.index("--start") + 1]) if "--start" in a else 0
    end = int(a[a.index("--end") + 1]) if "--end" in a else None
    outpath = a[a.index("--out") + 1] if "--out" in a else OUT

    fai = load_fai(FAI)
    chr_map = load_chr_map(FA)
    rd = FaiReader(FA, fai, chr_map)
    rows = []
    with open(TSV) as fh:
        for i, line in enumerate(fh):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if i == 0 and f[1] == "start":
                continue
            rows.append(f)
    if end is None:
        end = len(rows)
    rows = rows[start:end]
    print(f"processing {len(rows)} hits [{start}:{end}]", flush=True)
    res = []
    n = 0
    for f in rows:
        chrom, s, e = f[0], int(f[1]), int(f[2])
        try:
            seq = rd.get(chrom, s, e)
        except KeyError:
            res.append(f + ["NO_SEQ", "NA", "NA", "NA", "NA", "NA", "NA"])
            continue
        cat, smi, imi, cmi, sma, ima, cma = classify(seq)
        res.append(f + [cat, f"{smi:.0f}", f"{imi:.3f}", f"{cmi:.3f}",
                        f"{sma:.0f}", f"{ima:.3f}", f"{cma:.3f}"])
        n += 1
        if n % 500 == 0:
            print(f"  {n}/{len(rows)}", flush=True)
    with open(outpath, "w") as fh:
        fh.write("chrom\tstart\tend\tperiod_size\tcopies_aligned\tconsensus_size\t"
                 "match_percent\talignment_score\tconsensus_sequence\t"
                 "class\tsmi\tyimi\tcmi\tsma\tyima\tcma\n")
        for r in res:
            fh.write("\t".join(r) + "\n")
    print(f"wrote {outpath}", flush=True)


if __name__ == "__main__":
    main()
