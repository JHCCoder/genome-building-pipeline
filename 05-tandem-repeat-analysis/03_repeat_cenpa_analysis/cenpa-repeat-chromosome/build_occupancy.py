#!/usr/bin/env python
"""
build_occupancy.py — repeat-class occupancy around each chromosome's fixed
100-kb CENP-A core.

For every chromosome, orient relative to the 100-kb core (peak_window .. +100kb):
  short side = segment from core to NEARER chromosome end
  long  side = segment from core to FARTHER  chromosome end

Spatial intervals (mutually exclusive), all distances from the edge of the core:
  short distal       : 2–5 Mb
  short intermediate : 500 kb–2 Mb
  short proximal     : 0–500 kb
  core               : fixed 100-kb window
  long proximal      : 0–500 kb
  long intermediate  : 500 kb–2 Mb
  long distal        : 2–5 Mb
  background         : chromosome minus core minus all six flanking intervals

Intervals are clipped to chromosome boundaries; if an arm is shorter than a
requested interval, the available sequence is retained.

Repeat features: the 9 predefined TRF period bins, as pre-existing MERGED arrays
(period-enrichment/data/merged/arrays/bin{1..9}_arrays.bed; bedtools merge -d 0
per bin). Merging means overlapping annotations within a bin are already unique.

Per chromosome × interval × bin:
  interval_bp           actual bp after boundary clipping
  repeat_overlap_bp     unique bp of the bin overlapping the interval
  repeat_percent        = 100 * repeat_overlap_bp / interval_bp
  repeat_interval_count # of (merged) array intervals overlapping

Plus a QC table (partition check) and a genome-wide aggregate (chromosome="genome",
summed bp, not averaged).

Usage: python build_occupancy.py <work_dir> <out_prefix>
"""
import sys, os
import numpy as np
import pandas as pd

WORK = sys.argv[1] if len(sys.argv) > 1 else \
    "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
OUTP = sys.argv[2] if len(sys.argv) > 2 else os.path.join(WORK, "results", "repeat_occupancy")
# Optional 3rd arg: domain CSV with core_start/core_end columns (e.g. peak-based).
# Default: derive fixed 100-kb core from peak_window.
DOMAIN_CSV = sys.argv[3] if len(sys.argv) > 3 else None

SRC = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"
ARRAY_DIR = os.path.join(SRC, "period-enrichment", "data", "merged", "arrays")
DOMAINS = os.path.join(WORK, "data", "domains", "cenpa_domains_weighted.csv")
CHROMSIZES = os.path.join(SRC, "data", "chrom_sizes.txt")

CORE_WIN = 100000
BINS = list(range(1, 10))
BIN_LABEL = {1:"1-10 bp",2:"11-50 bp",3:"51-192 bp",4:"193-195 bp",5:"196-347 bp",
             6:"348-349 bp",7:"350-385 bp",8:"386-390 bp",9:"391+ bp"}

# ---- load chromosome sizes (chr1-28,X,Y) ----
sizes = {}
with open(CHROMSIZES) as f:
    for line in f:
        p = line.split()
        if p[0].startswith("chr"):
            sizes[p[0]] = int(p[1])
CHROMS = sorted(sizes.keys(), key=lambda c: (c[3:].isdigit(), int(c[3:]) if c[3:].isdigit() else (1 if c=="X" else 2)))

# ---- load domains ----
# If a domain CSV with explicit core_start/core_end is supplied (e.g. peak-based),
# use those coordinates directly. Otherwise derive a fixed 100-kb core anchored
# at peak_window.
dom = pd.read_csv(DOMAIN_CSV if DOMAIN_CSV else DOMAINS)
dom = dom[dom["chrom"].isin(CHROMS)].copy()
dom["chrom_len"] = dom["chrom"].map(sizes)
if DOMAIN_CSV and "core_start" in dom.columns and "core_end" in dom.columns:
    dom["core_start"] = dom["core_start"].astype(int)
    dom["core_end"]   = np.minimum(dom["core_end"].astype(int), dom["chrom_len"])
else:
    dom["core_start"] = dom["peak_window"]
    dom["core_end"]   = np.minimum(dom["peak_window"] + CORE_WIN, dom["chrom_len"])

# ---- load merged array BEDs per bin ----
def load_arrays(b):
    df = pd.read_csv(os.path.join(ARRAY_DIR, f"bin{b}_arrays.bed"), sep="\t",
                     header=None, names=["chrom","start","end","id","n_int"])
    return df[df["chrom"].isin(CHROMS)]
arrays = {b: load_arrays(b) for b in BINS}

# ---- spatial intervals ----
# edges (bp from core edge): proximal 0-500k, intermediate 500k-2M, distal 2-5M
EDGES = [("proximal", 0, 500000), ("intermediate", 500000, 2000000), ("distal", 2000000, 5000000)]
REGION_ORDER = {
    "short_distal":1, "short_intermediate":2, "short_proximal":3,
    "core":4,
    "long_proximal":5, "long_intermediate":6, "long_distal":7,
    "background":8,
}

def flank_interval_left(cs, lo, hi):
    """interval [lo,hi) bp from the core's LEFT edge, going left (toward 0)."""
    s = max(0, cs - hi)
    e = max(0, cs - lo)
    return (s, e) if e > s else None

def flank_interval_right(ce, clen, lo, hi):
    """interval [lo,hi) bp from the core's RIGHT edge, going right."""
    s = min(clen, ce + lo)
    e = min(clen, ce + hi)
    return (s, e) if e > s else None

def build_intervals(cs, ce, clen):
    """Return list of (side, region, istart, iend) for one chromosome.
    side: 'short'|'long'|'core'|'background'; region as labelled.
    All flank intervals are clipped to [0, clen]; empty intervals are dropped,
    so the seven focal regions + background partition the chromosome exactly."""
    left_len  = cs          # bp from core start to chromosome start (0)
    right_len = clen - ce   # bp from core end to chromosome end
    short_side = "left" if left_len <= right_len else "right"
    long_side  = "right" if short_side == "left" else "left"

    ivs = []
    for (region, lo, hi) in EDGES:
        # short side: nearest end
        if short_side == "left":
            iv = flank_interval_left(cs, lo, hi)
            if iv: ivs.append(("short", region, *iv))
        else:
            iv = flank_interval_right(ce, clen, lo, hi)
            if iv: ivs.append(("short", region, *iv))
        # long side: farthest end
        if long_side == "left":
            iv = flank_interval_left(cs, lo, hi)
            if iv: ivs.append(("long", region, *iv))
        else:
            iv = flank_interval_right(ce, clen, lo, hi)
            if iv: ivs.append(("long", region, *iv))

    # core
    ivs.append(("core", "core", cs, ce))

    # background = chromosome minus core minus all six flanks
    excl = [(s, e) for (_,_,s,e) in ivs]
    excl.sort()
    bg = []
    cur = 0
    for (s, e) in excl:
        if s > cur:
            bg.append((cur, s))
        cur = max(cur, e)
    if cur < clen:
        bg.append((cur, clen))
    for (s, e) in bg:
        if e > s:
            ivs.append(("background", "background", s, e))
    return ivs, short_side, long_side

# ---- overlap computation (numpy-based, vectorized) ----
def overlap_bp_and_count(arr, chrom, s, e):
    """unique bp overlap + number of array intervals overlapping [s,e) on chrom."""
    a = arr[arr["chrom"] == chrom]
    if len(a) == 0:
        return 0, 0
    st = a["start"].to_numpy(); en = a["end"].to_numpy()
    ov = np.maximum(0, np.minimum(en, e) - np.maximum(st, s))
    m = ov > 0
    return int(ov[m].sum()), int(m.sum())

rows = []
qc_rows = []
for _, d in dom.iterrows():
    chrom = d["chrom"]; cs = int(d["core_start"]); ce = int(d["core_end"])
    clen = int(d["chrom_len"])
    ivs, short_side, long_side = build_intervals(cs, ce, clen)

    # QC: interval lengths + partition check.
    # Note: background (and sometimes a flank) can split into multiple
    # non-contiguous segments; sum ALL segments for the partition check and
    # report the per-region bp as the SUM over its segments.
    valid = [(side, region, s, e) for (side, region, s, e) in ivs if e > s]
    total = sum(e - s for (_, _, s, e) in valid)
    reg_bp = {}
    for (side, region, s, e) in valid:
        key = f"{side}_{region}" if side not in ("core","background") else side
        reg_bp[key] = reg_bp.get(key, 0) + (e - s)
    qc_rows.append({
        "chromosome": chrom, "chrom_len": clen,
        "core_start": cs, "core_end": ce,
        "short_side": short_side, "long_side": long_side,
        "short_arm_bp": min(cs, clen - ce), "long_arm_bp": max(cs, clen - ce),
        "n_segments": len(valid),
        "partition_sum_bp": total,
        "partition_ok": total == clen,
        **reg_bp,
    })

    for (side, region, s, e) in ivs:
        ibp = e - s
        if ibp <= 0:
            continue
        region_key = f"{side}_{region}" if side not in ("core","background") else side
        for b in BINS:
            ov_bp, ov_n = overlap_bp_and_count(arrays[b], chrom, s, e)
            rows.append({
                "chromosome": chrom,
                "core_start": cs, "core_end": ce,
                "side": side, "region": region, "region_order": REGION_ORDER[region_key],
                "interval_start": s, "interval_end": e, "interval_bp": ibp,
                "repeat_bin": b, "repeat_bin_label": BIN_LABEL[b],
                "repeat_overlap_bp": ov_bp,
                "repeat_percent": 100.0 * ov_bp / ibp if ibp else 0.0,
                "repeat_interval_count": ov_n,
            })

occ = pd.DataFrame(rows)
qc = pd.DataFrame(qc_rows)

# ---- genome-wide aggregate (summed bp) ----
# NOTE: chrY is excluded from the genome-wide aggregate — no tandem-repeat (TRF)
# analysis was performed on chrY, so it contributes 0 repeat bp and would only
# dilute every genome-wide percentage. chrY remains in the per-chromosome rows.
occ_genome = (occ[occ["chromosome"] != "chrY"]
    .groupby(["side", "region", "region_order", "repeat_bin", "repeat_bin_label"], as_index=False)
    .agg(interval_bp=("interval_bp","sum"),
         repeat_overlap_bp=("repeat_overlap_bp","sum"),
         repeat_interval_count=("repeat_interval_count","sum")))
occ_genome["repeat_percent"] = 100.0 * occ_genome["repeat_overlap_bp"] / occ_genome["interval_bp"]
occ_genome["chromosome"] = "genome"
occ_genome.insert(0, "chromosome", occ_genome.pop("chromosome"))
occ_full = pd.concat([occ, occ_genome], ignore_index=True)

# ---- write outputs ----
occ.to_csv(OUTP + ".csv", index=False)
occ_full.to_csv(OUTP + "_with_genome.csv", index=False)
qc.to_csv(OUTP + "_qc.csv", index=False)

print(f"rows (per-chromosome): {len(occ)}")
print(f"rows (with genome aggregate): {len(occ_full)}")
print(f"QC rows: {len(qc)}; all partition_ok: {bool(qc['partition_ok'].all())}")
print("\n--- genome-wide summary (percent) ---")
g = occ_genome.pivot_table(index=["side","region"], columns="repeat_bin_label",
                           values="repeat_percent", aggfunc="first")
print(g.round(2).to_string())
print(f"\nwrote: {OUTP}.csv, {OUTP}_with_genome.csv, {OUTP}_qc.csv")
