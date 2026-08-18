#!/usr/bin/env python
"""
1a_compute_nh_weighted.py -- compute NH (alignments per read pair) and emit
1/NH-weighted fragment density in 100 kb windows.

Rationale: the k=100 BAMs carry no NH tags. NH = number of reported
alignments for the fragment. We count alignment records per read name
over the whole BAM (all secondary alignments, proper pairs only), divide
by 2 (both mates of a pair share the QNAME), then in a second pass weight
each fragment placement by 1/NH and accumulate it into the 100 kb window
containing the fragment midpoint.

Usage:
  python 1a_compute_nh_weighted.py <sample> <bam> <chrom_sizes> <outdir> <threads>

Outputs (in <outdir>):
  <sample>_nh.tsv.gz               qname -> NH (fragment-level, NH>=1)
  <sample>_win100kb_weighted.tsv   chrom  start  end  weighted_count
  <sample>_nh_hist.tsv             NH -> n_fragments (QC)
"""
import sys, gzip, collections
import pysam

SAMPLE, BAM, CHROM_SIZES, OUTDIR, THREADS = sys.argv[1:6]
THREADS = int(THREADS)

CHROMS = set()
with open(CHROM_SIZES) as f:
    for line in f:
        c = line.split()[0]
        if c.startswith("chr"):
            CHROMS.add(c)

def proper_pair(rec):
    # mapped, proper pair, not supplementary; include secondary (0x100)
    return (not rec.is_unmapped) and rec.is_proper_pair and (not rec.is_supplementary)

# ---------------------------------------------------------------------------
# Pass 1: count alignment records per read name
# ---------------------------------------------------------------------------
nh = collections.Counter()
n_records = 0
bam = pysam.AlignmentFile(BAM, "rb", threads=THREADS)
for rec in bam.fetch(until_eof=True):
    if not proper_pair(rec):
        continue
    nh[rec.query_name] += 1
    n_records += 1
bam.close()
print(f"[{SAMPLE}] pass1 records: {n_records}, distinct qnames: {len(nh)}", flush=True)

# Convert per-read counts to fragment NH (both mates share QNAME).
# Guard: if a QNAME has an odd count (e.g. a singleton mate), NH = ceil(count/2).
nh_frag = {q: max(1, (c + 1) // 2) for q, c in nh.items()}
del nh

# Write NH table + histogram
with gzip.open(f"{OUTDIR}/{SAMPLE}_nh.tsv.gz", "wt") as f:
    for q, c in nh_frag.items():
        f.write(f"{q}\t{c}\n")
hist = collections.Counter(nh_frag.values())
with open(f"{OUTDIR}/{SAMPLE}_nh_hist.tsv", "w") as f:
    f.write("NH\tn_fragments\n")
    for k in sorted(hist):
        f.write(f"{k}\t{hist[k]}\n")

# ---------------------------------------------------------------------------
# Pass 2: accumulate 1/NH-weighted fragment midpoints into 100 kb windows
# ---------------------------------------------------------------------------
win = collections.defaultdict(float)
W = 100000
n_placed = 0
bam = pysam.AlignmentFile(BAM, "rb", threads=THREADS)
for rec in bam.fetch(until_eof=True):
    if not proper_pair(rec):
        continue
    # Count each fragment placement once: use the first-in-pair (read1) records.
    if not (rec.is_read1):
        continue
    q = rec.query_name
    w = 1.0 / nh_frag[q]
    # fragment span from TLEN (ISIZE): midpoint = start + |TLEN|/2
    tlen = abs(rec.template_length) if rec.template_length is not None else 0
    if tlen > 0 and tlen < 1_000_000:  # sane fragment length
        mid = rec.reference_start + tlen // 2
    else:
        mid = rec.reference_start + (rec.reference_length or 0) // 2
    if rec.reference_name not in CHROMS:
        continue
    win[(rec.reference_name, mid // W * W)] += w
    n_placed += 1
bam.close()

with open(f"{OUTDIR}/{SAMPLE}_win100kb_weighted.tsv", "w") as f:
    f.write("chrom\tstart\tend\tweighted_count\n")
    for (chrom, s), v in sorted(win.items()):
        f.write(f"{chrom}\t{s}\t{s + W}\t{v:.4f}\n")

print(f"[{SAMPLE}] pass2 placed: {n_placed} fragment placements", flush=True)
print(f"[{SAMPLE}] total weighted: {sum(win.values()):.1f}", flush=True)
print(f"[{SAMPLE}] DONE", flush=True)
