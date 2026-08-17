#!/usr/bin/env python
"""make_occupancy_notebook.py — assemble the repeat-class occupancy notebook."""
import nbformat as nbf
import json, os

WORK = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
OUT = os.path.join(WORK, "repeat_occupancy_around_cenpa_core.ipynb")

nb = nbf.v4.new_notebook()
nb["metadata"] = {
    "kernelspec": {
        "display_name": "Python (python-visualizations)",
        "language": "python",
        "name": "python3",
    },
    "language_info": {"name": "python", "version": "3.13.1"},
}

cells = []

# ------------------------------------------------------------------ markdown
cells.append(nbf.v4.new_markdown_cell(
"""# Repeat-class occupancy around each chromosome's 100-kb CENP-A core

**Question.** For each degu chromosome, orient the analysis relative to the fixed
100-kb CENP-A core (anchored at the per-chromosome CENP-A peak window) and measure,
in mutually exclusive spatial intervals on the short and long side of the core,
how much of each of the **9 predefined TRF period bins** occupies the sequence.

**Repeat features.** The 9 tandem-repeat period bins, using the pre-existing
**merged arrays** (`period-enrichment/data/merged/arrays/bin{1..9}_arrays.bed`;
`bedtools merge -d 0` per bin). Merging means overlapping annotations within a bin
are counted once (no double-counting of bp).

**Spatial intervals** (distances from the edge of the CENP-A core; all clipped to
chromosome boundaries; if an arm is shorter than a requested interval the available
sequence is retained):

| side | region | distance from core edge |
|---|---|---|
| short | distal | 2–5 Mb |
| short | intermediate | 500 kb–2 Mb |
| short | proximal | 0–500 kb |
| — | **core** | fixed 100-kb window |
| long | proximal | 0–500 kb |
| long | intermediate | 500 kb–2 Mb |
| long | distal | 2–5 Mb |
| — | background | chromosome minus core minus all six flanks |

**Per chromosome × interval × bin:** `interval_bp`, `repeat_overlap_bp`
(unique bp, merged within bin), `repeat_percent` = 100·overlap/interval,
`repeat_interval_count` (# merged array intervals overlapping).

**Outputs**
- `results/repeat_occupancy.csv` — per-chromosome tidy dataframe
- `results/repeat_occupancy_with_genome.csv` — plus genome-wide aggregate rows (`chromosome = "genome"`, summed bp, not averaged)
- `results/repeat_occupancy_qc.csv` — one row per chromosome: length, core coords, short/long arm, and bp per interval; verifies the 7 focal regions + background partition each chromosome exactly
"""))

# ------------------------------------------------------------------ setup
cells.append(nbf.v4.new_code_cell(
"""import os
import numpy as np
import pandas as pd

WORK = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
SRC  = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment"

ARRAY_DIR   = os.path.join(SRC, "period-enrichment", "data", "merged", "arrays")
DOMAINS     = os.path.join(WORK, "data", "domains", "cenpa_domains_weighted.csv")
CHROMSIZES  = os.path.join(SRC, "data", "chrom_sizes.txt")
OUT_PREFIX  = os.path.join(WORK, "results", "repeat_occupancy")

CORE_WIN = 100_000
BINS     = list(range(1, 10))
BIN_LABEL = {1:"1-10 bp",2:"11-50 bp",3:"51-192 bp",4:"193-195 bp",5:"196-347 bp",
             6:"348-349 bp",7:"350-385 bp",8:"386-390 bp",9:"391+ bp"}

# ---- chromosome sizes (chr1-28, X, Y) ----
sizes = {}
with open(CHROMSIZES) as f:
    for line in f:
        p = line.split()
        if p[0].startswith("chr"):
            sizes[p[0]] = int(p[1])
CHROMS = sorted(sizes, key=lambda c: (c[3:].isdigit(), int(c[3:]) if c[3:].isdigit() else (1 if c=="X" else 2)))
print(f"{len(CHROMS)} chromosomes, total {sum(sizes.values())/1e6:.0f} Mb")
"""))

# ------------------------------------------------------------------ load domains
cells.append(nbf.v4.new_code_cell(
"""# ---- CENP-A domains; fixed 100-kb core anchored at the per-chromosome peak window ----
dom = pd.read_csv(DOMAINS)
dom = dom[dom["chrom"].isin(CHROMS)].copy()
dom["chrom_len"]  = dom["chrom"].map(sizes)
dom["core_start"] = dom["peak_window"]
dom["core_end"]   = np.minimum(dom["peak_window"] + CORE_WIN, dom["chrom_len"])
dom[["chrom","peak_window","core_start","core_end","chrom_len"]].head(8)
"""))

# ------------------------------------------------------------------ load arrays
cells.append(nbf.v4.new_code_cell(
"""# ---- load pre-existing MERGED array BEDs (one per period bin) ----
def load_arrays(b):
    df = pd.read_csv(os.path.join(ARRAY_DIR, f"bin{b}_arrays.bed"), sep="\\t",
                     header=None, names=["chrom","start","end","id","n_int"])
    return df[df["chrom"].isin(CHROMS)]

arrays = {b: load_arrays(b) for b in BINS}
for b in BINS:
    print(f"bin{b} ({BIN_LABEL[b]:>12s}): {len(arrays[b]):>7,} merged arrays")
"""))

# ------------------------------------------------------------------ intervals
cells.append(nbf.v4.new_code_cell(
"""# ---- spatial interval construction ----
# edges (bp from core edge): proximal 0-500k, intermediate 500k-2M, distal 2-5M
EDGES = [("proximal", 0, 500_000),
         ("intermediate", 500_000, 2_000_000),
         ("distal", 2_000_000, 5_000_000)]
REGION_ORDER = {
    "short_distal":1, "short_intermediate":2, "short_proximal":3,
    "core":4,
    "long_proximal":5, "long_intermediate":6, "long_distal":7,
    "background":8,
}

def flank_left(cs, lo, hi):
    s = max(0, cs - hi); e = max(0, cs - lo)
    return (s, e) if e > s else None

def flank_right(ce, clen, lo, hi):
    s = min(clen, ce + lo); e = min(clen, ce + hi)
    return (s, e) if e > s else None

def build_intervals(cs, ce, clen):
    \"\"\"Mutually-exclusive intervals for one chromosome.
    Returns (intervals, short_side, long_side).\"\"\"
    left_len  = cs
    right_len = clen - ce
    short_side = "left" if left_len <= right_len else "right"
    long_side  = "right" if short_side == "left" else "left"

    ivs = []
    for (region, lo, hi) in EDGES:
        for side in ("short", "long"):
            arm = short_side if side == "short" else long_side
            iv = (flank_left(cs, lo, hi) if arm == "left"
                  else flank_right(ce, clen, lo, hi))
            if iv: ivs.append((side, region, *iv))
    ivs.append(("core", "core", cs, ce))

    # background = chromosome minus core minus all six flanks
    excl = sorted((s, e) for (_, _, s, e) in ivs)
    bg, cur = [], 0
    for (s, e) in excl:
        if s > cur: bg.append((cur, s))
        cur = max(cur, e)
    if cur < clen: bg.append((cur, clen))
    for (s, e) in bg:
        if e > s: ivs.append(("background", "background", s, e))
    return ivs, short_side, long_side
"""))

# ------------------------------------------------------------------ compute
cells.append(nbf.v4.new_code_cell(
"""# ---- vectorized overlap: unique bp + # intervals overlapping [s,e) ----
def overlap_bp_and_count(arr, chrom, s, e):
    a = arr[arr["chrom"] == chrom]
    if len(a) == 0:
        return 0, 0
    st = a["start"].to_numpy(); en = a["end"].to_numpy()
    ov = np.maximum(0, np.minimum(en, e) - np.maximum(st, s))
    m = ov > 0
    return int(ov[m].sum()), int(m.sum())

rows, qc_rows = [], []
for _, d in dom.iterrows():
    chrom = d["chrom"]; cs = int(d["core_start"]); ce = int(d["core_end"])
    clen  = int(d["chrom_len"])
    ivs, short_side, long_side = build_intervals(cs, ce, clen)

    # ---- QC row (sum over ALL segments; a region may be multi-segment) ----
    valid = [(side, reg, s, e) for (side, reg, s, e) in ivs if e > s]
    total = sum(e - s for (_, _, s, e) in valid)
    reg_bp = {}
    for (side, reg, s, e) in valid:
        key = f"{side}_{reg}" if side not in ("core","background") else side
        reg_bp[key] = reg_bp.get(key, 0) + (e - s)
    qc_rows.append({
        "chromosome": chrom, "chrom_len": clen,
        "core_start": cs, "core_end": ce,
        "short_side": short_side, "long_side": long_side,
        "short_arm_bp": min(cs, clen - ce), "long_arm_bp": max(cs, clen - ce),
        "n_segments": len(valid), "partition_sum_bp": total,
        "partition_ok": total == clen, **reg_bp,
    })

    # ---- occupancy rows ----
    for (side, region, s, e) in ivs:
        ibp = e - s
        if ibp <= 0: continue
        region_key = f"{side}_{region}" if side not in ("core","background") else side
        for b in BINS:
            ov_bp, ov_n = overlap_bp_and_count(arrays[b], chrom, s, e)
            rows.append({
                "chromosome": chrom, "core_start": cs, "core_end": ce,
                "side": side, "region": region, "region_order": REGION_ORDER[region_key],
                "interval_start": s, "interval_end": e, "interval_bp": ibp,
                "repeat_bin": b, "repeat_bin_label": BIN_LABEL[b],
                "repeat_overlap_bp": ov_bp,
                "repeat_percent": 100.0 * ov_bp / ibp if ibp else 0.0,
                "repeat_interval_count": ov_n,
            })

occ = pd.DataFrame(rows)
qc  = pd.DataFrame(qc_rows)
print(f"per-chromosome rows: {len(occ)}")
print(f"QC rows: {len(qc)}; all partition_ok = {qc['partition_ok'].all()}")
"""))

# ------------------------------------------------------------------ genome aggregate
cells.append(nbf.v4.new_code_cell(
"""# ---- genome-wide aggregate: sum bp across chromosomes (NOT averaged) ----
occ_genome = (occ
    .groupby(["side","region","region_order","repeat_bin","repeat_bin_label"], as_index=False)
    .agg(interval_bp=("interval_bp","sum"),
         repeat_overlap_bp=("repeat_overlap_bp","sum"),
         repeat_interval_count=("repeat_interval_count","sum")))
occ_genome["repeat_percent"] = 100.0 * occ_genome["repeat_overlap_bp"] / occ_genome["interval_bp"]
occ_genome.insert(0, "chromosome", "genome")
occ_full = pd.concat([occ, occ_genome], ignore_index=True)
print(f"with genome rows: {len(occ_full)}")
"""))

# ------------------------------------------------------------------ save
cells.append(nbf.v4.new_code_cell(
"""# ---- write outputs ----
occ.to_csv(OUT_PREFIX + ".csv", index=False)
occ_full.to_csv(OUT_PREFIX + "_with_genome.csv", index=False)
qc.to_csv(OUT_PREFIX + "_qc.csv", index=False)
print("wrote:")
print("  ", OUT_PREFIX + ".csv")
print("  ", OUT_PREFIX + "_with_genome.csv")
print("  ", OUT_PREFIX + "_qc.csv")
"""))

# ------------------------------------------------------------------ QC view
cells.append(nbf.v4.new_code_cell(
"""# ---- QC: partition check ----
qc[["chromosome","chrom_len","core_start","core_end","short_side","long_side",
    "short_arm_bp","long_arm_bp","n_segments","partition_sum_bp","partition_ok"]]
"""))

# ------------------------------------------------------------------ genome summary
cells.append(nbf.v4.new_code_cell(
"""# ---- genome-wide 9-bin composition by spatial region (% of interval bp) ----
g = occ_genome.pivot_table(index=["side","region"], columns="repeat_bin_label",
                           values="repeat_percent", aggfunc="first")
g.round(2)
"""))

# ------------------------------------------------------------------ per-chr table
cells.append(nbf.v4.new_code_cell(
"""# ---- per-chromosome 349-bp (bin6) occupancy by region ----
p6 = occ[occ["repeat_bin"] == 6].pivot_table(
        index="chromosome", columns=["side","region"], values="repeat_percent", aggfunc="first")
# reorder regions short-distal ... background
order = [(s, r) for s, r in
         [("short","distal"),("short","intermediate"),("short","proximal"),
          ("core","core"),("long","proximal"),("long","intermediate"),
          ("long","distal"),("background","background")]]
p6 = p6[[c for c in order if c in p6.columns]].round(1)
p6.head(30)
"""))

# ------------------------------------------------------------------ signal check
cells.append(nbf.v4.new_code_cell(
"""# ---- sanity: chr4 core must be 100% 349-bp (peak sits inside the 8.7 Mb 349 array) ----
occ[(occ["chromosome"]=="chr4") & (occ["region"]=="core") & (occ["repeat_bin"]==6)] \\
    [["interval_bp","repeat_overlap_bp","repeat_percent","repeat_interval_count"]]
"""))

nb["cells"] = cells
with open(OUT, "w") as f:
    nbf.write(nb, f)
print("wrote", OUT)
