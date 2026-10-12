#!/usr/bin/env python3
"""
Figure 4D/E, step 1: genome-wide signal per chromosome-aligned interval around the
dominant CENP-A-enriched region (peak-based core), written to one tidy cache table.

Same calculations as repeat_occupancy_genomewide_plot.ipynb (cells 2-3, 8-11, 13-14, 16):
  * intervals per chromosome: the core (MACS2 CENP-A peak), then 0-500 kb, 0.5-2 Mb and
    2-5 Mb on the shorter and the longer arm, and the remaining background;
  * CENP-A  = mean of XG_150/XG_151 1/NH-weighted 100-kb window coverage;
    H3K27ac = XG_152 1/NH-weighted 100-kb window coverage;
  * HiFi    = mean 15-kb window depth, all reads (unfiltered) and primary alignments only;
  * TRF     = % of interval bp covered by merged arrays of a period bin (from the
              occupancy table; segments of the same interval are summed first).
Per chromosome, window values are overlap-weighted means within each interval; the
genome-wide value is the bp-weighted mean over chromosomes (chrY excluded: no TRF
analysis was run on chrY).

Run with the python-visualizations env:
    python 11_prepare_figure4DE_intervals.py
Writes $FIG4_INTERVALS_TSV (metric, x, interval, value).
"""

import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd


def _cfg(name):
    val = os.environ.get(name)
    if val:
        return val
    root = os.path.abspath(os.path.dirname(__file__))
    # repo-root config.sh (next to .git): this tree also has a local config.sh
    while not (os.path.isfile(os.path.join(root, "config.sh")) and os.path.exists(os.path.join(root, ".git"))) \
            and root != os.path.dirname(root):
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


WORK = _cfg("FIG4_INTERVAL_DIR")
INTERVAL_ORDER = [("short", "distal"), ("short", "intermediate"), ("short", "proximal"),
                  ("core", "core"),
                  ("long", "proximal"), ("long", "intermediate"), ("long", "distal"),
                  ("background", "background")]
INTERVAL_NAME = ["short 2-5 Mb", "short 0.5-2 Mb", "short 0-0.5 Mb", "core",
                 "long 0-0.5 Mb", "long 0.5-2 Mb", "long 2-5 Mb", "background"]
EDGES = [("proximal", 0, 500_000), ("intermediate", 500_000, 2_000_000), ("distal", 2_000_000, 5_000_000)]
TRF_BINS = ["1-10 bp", "193-195 bp", "348-349 bp", "386-390 bp"]

sizes = {}
with open(_cfg("FIG4_CHROM_SIZES")) as f:
    for line in f:
        p = line.split()
        if p[0].startswith("chr"):
            sizes[p[0]] = int(p[1])
CHROMS = [c for c in sizes if c != "chrY"]

dom = pd.read_csv(os.path.join(WORK, "data", "domains", "cenpa_domains_peakbased.csv"))
dom = dom[dom["chrom"].isin(CHROMS)].copy()
dom["chrom_len"] = dom["chrom"].map(sizes)
dom["core_start"] = dom["core_start"].astype(int)
dom["core_end"] = np.minimum(dom["core_end"].astype(int), dom["chrom_len"])


def build_intervals(cs, ce, clen):
    """(side, region, start, end) for the core, both flanks and the background."""
    left_len, right_len = cs, clen - ce
    short_side = "left" if left_len <= right_len else "right"
    ivs = []
    for region, lo, hi in EDGES:
        for side in ("short", "long"):
            arm = short_side if side == "short" else ("right" if short_side == "left" else "left")
            if arm == "left":
                s, e = max(0, cs - hi), max(0, cs - lo)
            else:
                s, e = min(clen, ce + lo), min(clen, ce + hi)
            if e > s:
                ivs.append((side, region, s, e))
    ivs.append(("core", "core", cs, ce))
    cur = 0
    for s, e in sorted((s, e) for (_, _, s, e) in ivs):
        if s > cur:
            ivs.append(("background", "background", cur, s))
        cur = max(cur, e)
    if cur < clen:
        ivs.append(("background", "background", cur, clen))
    return ivs


def interval_means(win, col):
    """Genome-wide bp-weighted mean of `col` per interval (win: chrom, start, end, col)."""
    rows = []
    for _, d in dom.iterrows():
        sub = win[win["chrom"] == d["chrom"]]
        for side, region, s, e in build_intervals(d["core_start"], d["core_end"], d["chrom_len"]):
            w = sub[(sub["start"] < e) & (sub["end"] > s)]
            if len(w) == 0:
                continue
            ov = np.maximum(0, np.minimum(w["end"], e) - np.maximum(w["start"], s))
            rows.append((side, region, e - s, (w[col] * ov).sum() / ov.sum()))
    t = pd.DataFrame(rows, columns=["side", "region", "bp", "mean"])
    g = t.groupby(["side", "region"]).apply(lambda x: (x["mean"] * x["bp"]).sum() / x["bp"].sum(),
                                           include_groups=False)
    return [g[k] for k in INTERVAL_ORDER]


out = []


def add(metric, values):
    for x, (name, v) in enumerate(zip(INTERVAL_NAME, values)):
        out.append({"metric": metric, "x": x, "interval": name, "value": v})


# CENP-A (mean of two replicates) and H3K27ac, 1/NH-weighted 100-kb windows
cov = None
for samp, col in [("XG_150", "s150"), ("XG_151", "s151"), ("XG_152", "h3")]:
    t = pd.read_csv(os.path.join(WORK, "data", "coverage_weighted", f"{samp}_win100kb_weighted.tsv"),
                    sep="\t", header=0, names=["chrom", "start", "end", col])
    cov = t if cov is None else cov.merge(t, on=["chrom", "start", "end"])
cov["cenpa"] = (cov["s150"] + cov["s151"]) / 2
cov = cov[cov["chrom"].isin(CHROMS)]
add("CENP-A", interval_means(cov, "cenpa"))
add("H3K27ac", interval_means(cov, "h3"))

# HiFi depth, 15-kb windows
for metric, var in [("HiFi all reads", "FIG4_HIFI_ALL_15KB"), ("HiFi primary", "FIG4_HIFI_PRIMARY_15KB")]:
    h = pd.read_csv(_cfg(var), sep="\t", header=None, names=["chrom", "start", "end", "depth"])
    add(metric, interval_means(h[h["chrom"].isin(CHROMS)], "depth"))

# TRF array occupancy (% of interval bp covered by merged arrays of the bin)
occ = pd.read_csv(os.path.join(WORK, "results", "repeat_occupancy_peakbased_with_genome.csv"))
gw = occ[occ["chromosome"] == "genome"]
for b in TRF_BINS:
    g = gw[gw["repeat_bin_label"] == b].groupby(["side", "region"])[["repeat_overlap_bp", "interval_bp"]].sum()
    add(b, [100 * g.loc[k, "repeat_overlap_bp"] / g.loc[k, "interval_bp"] for k in INTERVAL_ORDER])

res = pd.DataFrame(out)
path = _cfg("FIG4_INTERVALS_TSV")
os.makedirs(os.path.dirname(path), exist_ok=True)
if os.path.exists(path):
    os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
res.to_csv(path, sep="\t", index=False, float_format="%.6g")
print(res.pivot(index="x", columns="metric", values="value").round(3).to_string())
print("wrote", path)
