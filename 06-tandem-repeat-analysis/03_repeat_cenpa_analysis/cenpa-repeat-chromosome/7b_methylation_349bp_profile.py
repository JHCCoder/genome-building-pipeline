#!/usr/bin/env python
# ============================================================================
# 7b_methylation_349bp_profile.py -- CpG methylation of the 348-349 bp
# satellite, restricted to the satellite's own bases, in the SAME spatial
# intervals as the repeat-occupancy / all-CpG methylation analyses
# (6_methylation_profile.py, build_occupancy.py).
#
# Why 349-bp only: the 348-349 bp TRF period bin (bin 6) is the degu
# centromeric satellite. The all-CpG methylation panel (panel_methylation_
# around_cenpa_core) mixes satellite and non-satellite CpGs; this metric
# isolates the satellite so its methylation can be read off directly at each
# distance from the CENP-A core.
#
# Intervals (identical to 6_methylation_profile.py / build_occupancy.py):
#   core  = the per-chromosome peak-based CENP-A core (cenpa_domains_peakbased.csv)
#   short/long {distal 2-5 Mb, intermediate 500 kb-2 Mb, proximal 0-500 kb}
#   background = chromosome minus core minus all six flanks
#
# ARRAY SPLITTING: a merged 349-bp array frequently spans an interval boundary
# (e.g. chr1's 5-Mb satellite block crosses core/proximal/intermediate/distal).
# Every array is therefore CLIPPED at the interval edge -- each interval only
# receives the array bp / CpGs physically inside it -- so an array spanning N
# intervals contributes its proportional piece to each. This mirrors the
# overlap_bp_and_count() clip in build_occupancy.py, so bin6_overlap_bp here
# must exactly equal repeat_overlap_bp in repeat_occupancy_peakbased_with_genome.csv
# for bin 6 (verified at the end of this script).
#
# Per interval per chromosome (satellite CpGs only, WGBS cov >= MIN_COV):
#   n_cpgs_cov5      # covered satellite CpGs with cov >= MIN_COV
#   mean_frac        unweighted mean methylation over those CpGs
#   wmean_frac       coverage-weighted mean methylation (read-level proportion)
#   mean_cov         mean coverage over covered CpGs
#   bin6_overlap_bp  bp of 349-bp satellite inside the interval (after clip)
#   bin6_occupancy_pct = 100 * bin6_overlap_bp / interval_bp
# Intervals with < MIN_CPG covered CpGs are flagged sparse_ok = 0.
# Plus genome-wide pooled rows (weighted by n_cpgs across chromosomes with data).
#
# Usage: python 7b_methylation_349bp_profile.py <work_dir>   (python-visualizations env)
# ============================================================================
import os, sys
import numpy as np
import pandas as pd
import pyBigWig

WORK = sys.argv[1] if len(sys.argv) > 1 else \
    "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
DOM    = os.path.join(WORK, "data", "domains", "cenpa_domains_peakbased.csv")
BW_DIR = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/circos-plot/feature-overview"
FRAC   = os.path.join(BW_DIR, "degu_6834_PFC_1.CGN-both.frac.bw")
COV    = os.path.join(BW_DIR, "degu_6834_PFC_1.CGN-both.cov.bw")
B6     = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/period-enrichment/data/merged/arrays/bin6_arrays.bed"
OCC    = os.path.join(WORK, "results", "repeat_occupancy_peakbased_with_genome.csv")
OUT    = os.path.join(WORK, "results", "methylation_349bp_satellite_by_region.csv")

MIN_COV = 5      # per-CpG coverage filter
MIN_CPG = 10     # intervals with fewer covered satellite CpGs flagged as sparse

# ---- spatial intervals (mirrors 6_methylation_profile.py / build_occupancy.py) ----
EDGES = [("proximal", 0, 500_000),
         ("intermediate", 500_000, 2_000_000),
         ("distal", 2_000_000, 5_000_000)]
REGION_ORDER = {"short_distal": 1, "short_intermediate": 2, "short_proximal": 3,
                "core": 4,
                "long_proximal": 5, "long_intermediate": 6, "long_distal": 7,
                "background": 8}

def flank_left(cs, lo, hi):
    s = max(0, cs - hi); e = max(0, cs - lo)
    return (s, e) if e > s else None

def flank_right(ce, clen, lo, hi):
    s = min(clen, ce + lo); e = min(clen, ce + hi)
    return (s, e) if e > s else None

def build_intervals(cs, ce, clen):
    """Mutually-exclusive, exhaustive intervals for one chromosome (same as the
    occupancy / methylation notebooks): core + 3 flanks on each arm + background.
    Returns (side, region, start, end)."""
    left_len, right_len = cs, clen - ce
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
    excl = sorted((s, e) for (_, _, s, e) in ivs)
    bg, cur = [], 0
    for (s, e) in excl:
        if s > cur: bg.append((cur, s))
        cur = max(cur, e)
    if cur < clen: bg.append((cur, clen))
    for (s, e) in bg:
        if e > s: ivs.append(("background", "background", s, e))
    return ivs

# ---- load 349-bp (bin 6) merged arrays, sorted per chromosome ----
arr = {}
with open(B6) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 3: continue
        arr.setdefault(p[0], []).append((int(p[1]), int(p[2])))
for c in arr: arr[c].sort()

# ---- satellite methylation stats within one interval ----
# Each array is CLIPPED to the interval [s,e); only the clipped piece [a,b] is
# scanned for WGBS CpGs. This is the "split arrays at the interval boundary"
# rule: an array spanning several intervals contributes per-interval pieces.
def satellite_stats(fr, cv, asegs, chrom, s, e):
    fracs, covs, arrbp = [], [], 0
    for (as_, ae) in asegs:
        if as_ >= e: break                # sorted by start; no later array overlaps
        a, b = max(s, as_), min(e, ae)
        if b <= a: continue
        arrbp += b - a                    # clipped (split) bp inside this interval
        fi = fr.intervals(chrom, a, b) or []
        ci = cv.intervals(chrom, a, b) or []
        covd = {stt: v for (stt, _, v) in ci}
        for (stt, en, fv) in fi:
            c = covd.get(stt, np.nan)
            if not np.isnan(c) and c >= MIN_COV:
                fracs.append(fv); covs.append(c)
    n = len(fracs)
    return dict(n_cpgs_cov5=n,
                mean_frac=float(np.mean(fracs)) if n else np.nan,
                wmean_frac=float(np.sum(np.array(fracs) * np.array(covs)) / np.sum(covs)) if n else np.nan,
                mean_cov=float(np.mean(covs)) if n else np.nan,
                bin6_overlap_bp=arrbp,
                sparse_ok=1 if n >= MIN_CPG else 0)

def main():
    dom = pd.read_csv(DOM)
    fr = pyBigWig.open(FRAC); cv = pyBigWig.open(COV)

    rows = []
    for _, d in dom.iterrows():
        chrom = str(d["chrom"]); cs = int(d["core_start"]); ce = int(d["core_end"])
        clen = int(d["chrom_len"])
        asegs = arr.get(chrom, [])
        for (side, region, s, e) in build_intervals(cs, ce, clen):
            key = side if side in ("core", "background") else f"{side}_{region}"
            st = satellite_stats(fr, cv, asegs, chrom, s, e)
            ibp = e - s
            rows.append({
                "chromosome": chrom, "core_start": cs, "core_end": ce,
                "side": side, "region": region, "region_order": REGION_ORDER[key],
                "interval_start": s, "interval_end": e, "interval_bp": ibp,
                "bin6_occupancy_pct": 100 * st["bin6_overlap_bp"] / ibp if ibp else 0.0,
                **{k: st[k] for k in ("bin6_overlap_bp", "n_cpgs_cov5",
                                       "mean_frac", "wmean_frac", "mean_cov", "sparse_ok")}})
    fr.close(); cv.close()

    meth = pd.DataFrame(rows)

    # ---- validation: per chromosome, clipped 349-bp bp across the tiling
    # intervals must equal total 349-bp array bp on that chromosome ----
    tot_in = meth.groupby("chromosome")["bin6_overlap_bp"].sum()
    tot_arr = pd.Series({c: sum(ae - as_ for as_, ae in segs) for c, segs in arr.items()})
    bad = []
    for c in tot_arr.index:
        if c in tot_in and abs(tot_in[c] - tot_arr[c]) > 0:
            bad.append((c, tot_arr[c], tot_in.get(c, 0)))
    if bad:
        print("  !! PARTITION MISMATCH (chrom, total_array_bp, clipped_bp):")
        for r in bad: print("     ", r)
    else:
        print(f"  partition OK: clipped 349-bp bp tiles each chromosome exactly "
              f"({len(tot_arr)} chromosomes, no bp lost/doubled)")

    # ---- cross-check against the occupancy table (bin 6, same clip rule) ----
    # background (and sometimes a flank) splits into multiple segments, so
    # compare the per-(chrom, side, region) SUM of overlap bp between tables.
    occ = pd.read_csv(OCC)
    occ6 = occ[(occ["repeat_bin"] == 6)]
    mine = (meth.groupby(["chromosome", "side", "region"])["bin6_overlap_bp"].sum()
            .reset_index())
    occs = (occ6.groupby(["chromosome", "side", "region"])["repeat_overlap_bp"].sum()
            .reset_index())
    chk = mine.merge(occs, on=["chromosome", "side", "region"])
    diff = chk.loc[chk["bin6_overlap_bp"] != chk["repeat_overlap_bp"]]
    if len(diff):
        print(f"  !! bin6_overlap differs from occupancy table for {len(diff)} rows")
        print(diff.to_string(index=False))
    else:
        print(f"  cross-check OK: bin6_overlap_bp == repeat_occupancy_peakbased_with_genome.csv "
              f"(bin 6), all {len(chk)} (chrom, region) rows")

    # ---- genome-wide pooled rows (pooled across chromosomes WITH data,
    # weighted by n_cpgs_cov5 -- same convention as 6_methylation_profile.py) ----
    gw = []
    for (side, region, ro), g in meth.groupby(["side", "region", "region_order"]):
        allrows = meth[(meth["side"] == side) & (meth["region"] == region)]
        ibp = int(allrows["interval_bp"].sum()); obp = int(allrows["bin6_overlap_bp"].sum())
        gd = g[g["n_cpgs_cov5"] > 0]
        if len(gd) == 0:
            gw.append({"chromosome": "genome", "side": side, "region": region,
                       "region_order": ro, "interval_bp": ibp,
                       "bin6_overlap_bp": obp,
                       "bin6_occupancy_pct": 100 * obp / ibp if ibp else 0.0,
                       "n_cpgs_cov5": 0, "mean_frac": np.nan,
                       "wmean_frac": np.nan, "mean_cov": np.nan,
                       "sparse_ok": 0, "n_chr_with_cpg": 0})
            continue
        tot = int(gd["n_cpgs_cov5"].sum())
        gw.append({"chromosome": "genome", "side": side, "region": region,
                   "region_order": ro, "interval_bp": ibp,
                   "bin6_overlap_bp": obp,
                   "bin6_occupancy_pct": 100 * obp / ibp if ibp else 0.0,
                   "n_cpgs_cov5": tot,
                   "mean_frac": float(np.average(gd["mean_frac"], weights=gd["n_cpgs_cov5"])),
                   "wmean_frac": float(np.average(gd["wmean_frac"], weights=gd["n_cpgs_cov5"])),
                   "mean_cov": float(np.average(gd["mean_cov"], weights=gd["n_cpgs_cov5"])),
                   "sparse_ok": 1 if tot >= MIN_CPG else 0,
                   "n_chr_with_cpg": int((gd["n_cpgs_cov5"] > 0).sum())})
    meth = pd.concat([meth, pd.DataFrame(gw)], ignore_index=True)
    meth.to_csv(OUT, index=False)

    print(f"wrote {OUT}")
    view = meth[meth["chromosome"] == "genome"][
        ["region_order", "side", "region", "interval_bp", "bin6_overlap_bp",
         "bin6_occupancy_pct", "n_chr_with_cpg", "n_cpgs_cov5",
         "mean_frac", "wmean_frac", "mean_cov"]].sort_values("region_order")
    print(view.to_string(index=False))

if __name__ == "__main__":
    main()
