#!/usr/bin/env python
"""Figure 2B, step 1: per-window coverage of genes and repeats (cached TSV).

Same quantity as 02_genome_density_plot.ipynb: for each chromosome, the percentage of
each 1 Mb window covered by features (overlapping features merged first), then a
centred rolling mean over 10 windows. Features:
  Liftoff  - "gene" records of the Liftoff-transferred annotation
  BRAKER3  - "gene" records of the BRAKER3 de novo annotation
  Repeats  - all RepeatMasker records
Only chr1-chr28, chrX and chrY are kept.

Vectorised rewrite of the notebook's row-by-row loop (the repeat GFF is ~1 GB).
Intervals are treated as GFF 1-based closed [start, end] -> 0-based [start-1, end).

Run with the python-visualizations env (needs a few GB of RAM; not on a login node):
    python 02_prepare_figure2B_density.py
"""

import os
import subprocess
import sys

import numpy as np
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


WINDOW = 1_000_000
ROLL = 10                                   # windows in the centred rolling mean (notebook: 150000 // 15000)
CHROMS = [f"chr{i}" for i in range(1, 29)] + ["chrX", "chrY"]
TRACKS = [  # (track, config variable, feature type or None for all)
    ("Liftoff", "FIG2B_LIFTOFF_GFF", "gene"),
    ("BRAKER3", "FIG2B_BRAKER_GFF", "gene"),
    ("Repeats", "FIG2B_REPEAT_GFF", None),
]

fai = pd.read_csv(_cfg("FIG2B_FAI"), sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"])
lengths = dict(zip(fai["chrom"], fai["length"]))
missing = [c for c in CHROMS if c not in lengths]
assert not missing, f"not in .fai: {missing}"


def covered_upto(starts, ends, x):
    """Bases covered by merged, sorted intervals [starts, ends) in [0, x) for each x."""
    seg = ends - starts
    cum = np.concatenate([[0], np.cumsum(seg)])
    i = np.searchsorted(starts, x, side="right")          # intervals starting before x
    last = np.clip(i - 1, 0, None)
    partial = np.where(i > 0, np.clip(np.minimum(x, ends[last]) - starts[last], 0, None), 0)
    return cum[np.clip(i - 1, 0, None)] * (i > 0) + partial


def merge(starts, ends):
    order = np.argsort(starts, kind="stable")
    s, e = starts[order], ends[order]
    e_run = np.maximum.accumulate(e)
    new = np.concatenate([[True], s[1:] > e_run[:-1]])    # a gap before this interval
    grp = np.cumsum(new) - 1
    return s[new], np.maximum.reduceat(e_run, np.flatnonzero(new)) if len(s) else e[:0]


rows = []
for track, var, ftype in TRACKS:
    path = _cfg(var)
    print(f"[{track}] reading {path}", flush=True)
    gff = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 2, 3, 4],
                      names=["chrom", "type", "start", "end"], dtype={"chrom": str, "type": str})
    if ftype:
        gff = gff[gff["type"] == ftype]
    gff = gff[gff["chrom"].isin(CHROMS)]
    print(f"[{track}] {len(gff):,} features on {gff['chrom'].nunique()} chromosomes", flush=True)
    for chrom in CHROMS:
        sub = gff[gff["chrom"] == chrom]
        s, e = merge(sub["start"].to_numpy(np.int64) - 1, sub["end"].to_numpy(np.int64))
        L = lengths[chrom]
        bounds = np.append(np.arange(0, L, WINDOW), L)
        cov = np.diff(covered_upto(s, e, bounds)) if len(s) else np.zeros(len(bounds) - 1)
        pct = cov / np.diff(bounds) * 100
        smooth = pd.Series(pct).rolling(ROLL, min_periods=1, center=True).mean().to_numpy()
        rows.append(pd.DataFrame({"track": track, "chrom": chrom, "start": bounds[:-1], "end": bounds[1:],
                                  "covered_bp": cov.astype(np.int64), "pct": pct, "pct_smoothed": smooth}))

out = _cfg("FIG2B_DENSITY_TSV")
os.makedirs(os.path.dirname(out), exist_ok=True)
df = pd.concat(rows, ignore_index=True)
df.to_csv(out, sep="\t", index=False, float_format="%.4f")
print(f"wrote {out}: {len(df):,} rows")
print(df.groupby("track")["pct"].agg(["mean", "max"]).round(1))
