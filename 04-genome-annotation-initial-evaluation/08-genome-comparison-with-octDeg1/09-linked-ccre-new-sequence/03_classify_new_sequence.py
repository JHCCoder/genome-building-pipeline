#!/usr/bin/env python3
"""
Step 3: classify every peak–gene linked peak (union of the four OctDeg2.0 samples) as

  shared            >= 50% of the peak in OctDeg2.0 sequence aligned to OctDeg1.0 (step 1), or
                    an unaligned candidate with exactly one OctDeg1.0 BLAST locus (step 2)
  resolved_repeat   unaligned candidate whose sequence occurs at >= 2 OctDeg1.0 loci
                    (present in OctDeg1.0 as other copies; OctDeg2.0 resolves this copy)
  partial           unaligned candidate with BLAST hits, none covering >= 80% of the peak
  new_sequence      unaligned candidate with no BLAST hit at all (sequence absent from OctDeg1.0)

A BLAST locus = HSPs with >= 90% identity that together cover >= 80% of the peak, on one
OctDeg1.0 sequence within 5 kb of each other.

Outputs:
  $FIG3D_CCRE_DIR/03_linked_peaks_classified.tsv        one row per union peak
  $FIG3D_CCRE_DIR/03_linked_peaks_<class>.bed           peaks per class (for motif analysis)
  results/03_classification_summary.tsv                  counts per class, union and per sample

Run with the python-visualizations env:
    python 03_classify_new_sequence.py
"""

import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


def _cfg(name):
    val = os.environ.get(name)
    if val:
        return val
    root = os.path.abspath(os.path.dirname(__file__))
    while not os.path.isfile(os.path.join(root, "config.sh")) and root != os.path.dirname(root):
        root = os.path.dirname(root)
    # Source config.sh in a subshell so nested $VAR references are expanded.
    val = subprocess.check_output(
        ["bash", "-c", f'source "$1" 2>/dev/null; printf "%s" "${{{name}}}"', "_",
         os.path.join(root, "config.sh")],
        text=True,
    )
    if not val:
        sys.exit(f"ERROR: {name} is not set in the environment or in config.sh")
    return val


MIN_PID, MIN_QCOV, LOCUS_GAP = 90.0, 0.80, 5000
SAMPLES = ["181_PFC", "181_dHIP", "6997_dHIP", "7000_PFC"]
OUT = Path(_cfg("FIG3D_CCRE_DIR"))
RES = Path(__file__).resolve().parent / "results"
RES.mkdir(exist_ok=True)

# --- union peaks with aligned fraction (step 1) ----------------------------------------
union = pd.read_csv(OUT / "01_union.alignedFrac.bed", sep="\t", header=None)
union = union.iloc[:, [0, 1, 2, 3, 4, union.shape[1] - 1]]
union.columns = ["chrom", "start", "end", "n_samples", "samples", "aligned_frac"]
union["peak"] = union["chrom"] + ":" + union["start"].astype(str) + "-" + union["end"].astype(str)

# --- BLAST loci per candidate (step 2) --------------------------------------------------
hits = pd.read_csv(OUT / "02_candidates_vs_octdeg1.blastn.tsv", sep="\t", header=None,
                   names=["q", "s", "pid", "len", "qs", "qe", "ss", "se", "ev", "bits", "qlen"])
hits = hits[hits["pid"] >= MIN_PID]
n_any_hit = hits.groupby("q").size()


def n_loci(df):
    """Group HSPs into loci on each subject; count loci covering >= MIN_QCOV of the query."""
    qlen = df["qlen"].iat[0]
    loci = 0
    for _, g in df.groupby("s"):
        g = g.assign(lo=g[["ss", "se"]].min(axis=1)).sort_values("lo")
        clusters, cur, last = [], [], None
        for r in g.itertuples():
            if last is not None and r.lo - last > LOCUS_GAP:
                clusters.append(cur); cur = []
            cur.append(r); last = max(r.ss, r.se)
        clusters.append(cur)
        for c in clusters:
            covered = set()
            for r in c:
                covered.update(range(r.qs, r.qe + 1))
            if len(covered) / qlen >= MIN_QCOV:
                loci += 1
    return loci


full_loci = (hits.groupby("q")[["s", "pid", "len", "qs", "qe", "ss", "se", "qlen"]].apply(n_loci)
             if len(hits) else pd.Series(dtype=int))

def classify(r):
    if r.aligned_frac >= 0.5:
        return "shared"
    k = int(full_loci.get(r.peak, 0))
    if k == 1:
        return "shared"
    if k >= 2:
        return "resolved_repeat"
    return "partial" if r.peak in n_any_hit.index else "new_sequence"

union["octdeg1_loci"] = union["peak"].map(full_loci).fillna(0).astype(int)
union.loc[union["aligned_frac"] >= 0.5, "octdeg1_loci"] = pd.NA
union["class"] = union.apply(classify, axis=1)
union.to_csv(OUT / "03_linked_peaks_classified.tsv", sep="\t", index=False)
for c, g in union.groupby("class"):
    g[["chrom", "start", "end", "peak", "n_samples", "samples"]].to_csv(
        OUT / f"03_linked_peaks_{c}.bed", sep="\t", header=False, index=False)

# --- summary: union and per sample (a sample counts a union peak it contributed to) ------
order = ["shared", "resolved_repeat", "partial", "new_sequence"]
rows = [{"set": "union", **union["class"].value_counts().reindex(order, fill_value=0).to_dict()}]
for s in SAMPLES:
    sub = union[union["samples"].str.split(",").apply(lambda xs: s in xs)]
    rows.append({"set": s, **sub["class"].value_counts().reindex(order, fill_value=0).to_dict()})
summ = pd.DataFrame(rows)
summ["total"] = summ[order].sum(axis=1)
for c in order[1:]:
    summ[f"pct_{c}"] = (100 * summ[c] / summ["total"]).round(2)
summ.to_csv(RES / "03_classification_summary.tsv", sep="\t", index=False)
print(summ.to_string(index=False))
