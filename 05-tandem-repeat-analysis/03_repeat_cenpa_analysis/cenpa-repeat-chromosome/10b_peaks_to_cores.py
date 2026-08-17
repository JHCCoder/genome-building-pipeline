#!/usr/bin/env python
"""
10b_peaks_to_cores.py — redefine each chromosome's CENP-A core from MACS2 peaks.

For each chromosome, the new core is the MACS2 peak that overlaps the
strict-signal peak window (the chromosome's CENP-A anchor, from
cenpa_domains_weighted.csv peak_window). If several peaks overlap the anchor
window, the one with the highest signal (col 7 in narrowPeak) is used. If no
peak overlaps, fall back to the fixed 100 kb anchor window.

Usage: python 10b_peaks_to_cores.py <work_dir>
"""
import os, sys
import numpy as np
import pandas as pd

WORK = sys.argv[1] if len(sys.argv) > 1 else \
    "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
PEAKS = os.path.join(WORK, "data", "macs2_peaks", "cenpa_peaks.narrowPeak")
DOM   = os.path.join(WORK, "data", "domains", "cenpa_domains_weighted.csv")
OUT   = os.path.join(WORK, "data", "domains", "cenpa_domains_peakbased.csv")

CORE_WIN = 100_000

# chromosome lengths (chr1-28, X)
CHROMSIZES = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/data/chrom_sizes.txt"
chr_len = {}
with open(CHROMSIZES) as f:
    for line in f:
        p = line.split()
        if p[0].startswith("chr"):
            chr_len[p[0]] = int(p[1])

dom = pd.read_csv(DOM)
dom = dom[dom["chrom"] != "chrY"]                      # chrY excluded (no TRF)
dom["chrom_len"] = dom["chrom"].map(chr_len)
dom["anchor_start"] = dom["peak_window"]
dom["anchor_end"]   = np.minimum(dom["peak_window"] + CORE_WIN, dom["chrom_len"])

# narrowPeak: chrom start end name score strand signal pval qval peak
cols = ["chrom", "start", "end", "name", "score", "strand",
        "signal", "pval", "qval", "peak_offset"]
peaks = pd.read_csv(PEAKS, sep="\t", header=None, names=cols)
peaks = peaks[peaks["chrom"].isin(dom["chrom"])].copy()

rows = []
for _, d in dom.iterrows():
    chrom = d["chrom"]; a_s = int(d["anchor_start"]); a_e = int(d["anchor_end"])
    ov = peaks[(peaks["chrom"] == chrom) & (peaks["start"] < a_e) & (peaks["end"] > a_s)]
    if len(ov) == 0:
        # no peak overlaps the anchor -> fall back to the 100 kb anchor window
        rows.append({**d.to_dict(),
                     "core_start": a_s, "core_end": a_e, "core_size": a_e - a_s,
                     "peak_start": np.nan, "peak_end": np.nan, "peak_signal": np.nan,
                     "peak_source": "fallback_anchor"})
        continue
    best = ov.loc[ov["signal"].idxmax()]
    rows.append({**d.to_dict(),
                 "core_start": int(best["start"]), "core_end": int(best["end"]),
                 "core_size": int(best["end"] - best["start"]),
                 "peak_start": int(best["start"]), "peak_end": int(best["end"]),
                 "peak_signal": float(best["signal"]),
                 "peak_source": "macs2_peak"})

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print(f"wrote {OUT}")
print(out[["chrom", "anchor_start", "anchor_end",
           "core_start", "core_end", "core_size", "peak_source"]].to_string(index=False))
print(f"\nmean core size: {out['core_size'].mean()/1000:.1f} kb; "
      f"median: {out['core_size'].median()/1000:.1f} kb; "
      f"peak-based: {(out['peak_source']=='macs2_peak').sum()}/{len(out)} chromosomes")
