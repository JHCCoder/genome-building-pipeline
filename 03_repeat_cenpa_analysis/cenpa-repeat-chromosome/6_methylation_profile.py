#!/usr/bin/env python
# ============================================================================
# 6_methylation_profile.py -- CpG methylation level around the peak-based
# CENP-A cores, in the SAME spatial intervals as the repeat-occupancy analysis
# (repeat_occupancy_around_cenpa_core.ipynb).
#
# Intervals (identical construction to the occupancy notebook):
#   core  = the per-chromosome peak-based CENP-A core (cenpa_domains_peakbased.csv)
#   short/long {distal 2-5 Mb, intermediate 500 kb-2 Mb, proximal 0-500 kb}
#   background = chromosome minus core minus all six flanks
# short/long sides are the shorter/longer arms measured from the core edge.
#
# Methylation data (WGBS, CpG = CGN context, both strands):
#   degu_6834_PFC_1.CGN-both.frac.bw   -- fraction methylated (0..1) per CpG
#   degu_6834_PFC_1.CGN-both.cov.bw    -- read coverage per CpG (for filtering)
#
# Per interval per chromosome:
#   n_cpgs_cov5      # covered CpGs with cov >= MIN_COV
#   mean_frac        unweighted mean methylation over those CpGs
#   wmean_frac       coverage-weighted mean methylation (read-level proportion)
#   mean_cov         mean coverage over covered CpGs
# Intervals with < MIN_CPG covered CpGs are flagged (sparse_ok = 0).
#
# Usage: python 6_methylation_profile.py <work_dir>   (python-visualizations env)
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
OUT    = os.path.join(WORK, "results", "methylation_around_cenpa_core.csv")

MIN_COV = 5      # per-CpG coverage filter
MIN_CPG = 10     # intervals with fewer covered CpGs flagged as sparse

# ---- spatial intervals (mirrors repeat_occupancy_around_cenpa_core.ipynb) ----
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
    """Mutually-exclusive intervals for one chromosome (same as notebook)."""
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

# ---- per-interval methylation stats ----
def interval_stats(fr, cv, chrom, s, e):
    fi = fr.intervals(chrom, s, e) or []
    ci = cv.intervals(chrom, s, e) or []
    # build cov dict keyed by start (cov and frac bigwigs share CpG positions)
    covd = {st: v for (st, _, v) in ci}
    fracs, covs = [], []
    for (st, en, fv) in fi:
        c = covd.get(st, np.nan)
        if np.isnan(c) or c < MIN_COV:
            continue
        fracs.append(fv); covs.append(c)
    n = len(fracs)
    return dict(n_cpgs_cov5=n,
                mean_frac=float(np.mean(fracs)) if n else np.nan,
                wmean_frac=float(np.sum(np.array(fracs) * np.array(covs)) / np.sum(covs)) if n else np.nan,
                mean_cov=float(np.mean(covs)) if n else np.nan,
                sparse_ok=1 if n >= MIN_CPG else 0)

def main():
    dom = pd.read_csv(DOM)
    fr = pyBigWig.open(FRAC); cv = pyBigWig.open(COV)

    rows = []
    for _, d in dom.iterrows():
        chrom = str(d["chrom"]); cs = int(d["core_start"]); ce = int(d["core_end"])
        clen = int(d["chrom_len"])
        for (side, region, s, e) in build_intervals(cs, ce, clen):
            key = side if side in ("core", "background") else f"{side}_{region}"
            st = interval_stats(fr, cv, chrom, s, e)
            rows.append({
                "chromosome": chrom, "core_start": cs, "core_end": ce,
                "side": side, "region": region, "region_order": REGION_ORDER[key],
                "interval_start": s, "interval_end": e, "interval_bp": e - s,
                **st})
    fr.close(); cv.close()

    meth = pd.DataFrame(rows)
    meth.to_csv(OUT, index=False)

    # genome-wide pooled rows (sum bp / pooled CpGs across chromosomes, like the
    # occupancy notebook: pooled, not averaged)
    gw = []
    for (side, region, ro), g in meth.groupby(["side", "region", "region_order"]):
        g = g[g["n_cpgs_cov5"] > 0]
        if len(g) == 0:
            gw.append({"chromosome": "genome", "side": side, "region": region,
                       "region_order": ro, "interval_bp": int(g["interval_bp"].sum()),
                       "n_cpgs_cov5": 0, "mean_frac": np.nan, "wmean_frac": np.nan,
                       "mean_cov": np.nan, "sparse_ok": 0})
            continue
        tot_cpgs = int(g["n_cpgs_cov5"].sum())
        # pooled mean = weighted by n_cpgs (each CpG contributes once)
        mean_frac  = float(np.average(g["mean_frac"], weights=g["n_cpgs_cov5"]))
        wmean_frac = float(np.average(g["wmean_frac"], weights=g["n_cpgs_cov5"]))
        mean_cov   = float(np.average(g["mean_cov"], weights=g["n_cpgs_cov5"]))
        gw.append({"chromosome": "genome", "side": side, "region": region,
                   "region_order": ro, "interval_bp": int(g["interval_bp"].sum()),
                   "n_cpgs_cov5": tot_cpgs, "mean_frac": mean_frac,
                   "wmean_frac": wmean_frac, "mean_cov": mean_cov,
                   "sparse_ok": 1 if tot_cpgs >= MIN_CPG else 0})
    meth = pd.concat([meth, pd.DataFrame(gw)], ignore_index=True)
    meth.to_csv(OUT, index=False)

    print(f"wrote {OUT}  ({len(dom)} chromosomes x intervals + genome rows)")
    view = meth[meth["chromosome"] == "genome"][
        ["region_order", "side", "region", "interval_bp", "n_cpgs_cov5",
         "mean_frac", "wmean_frac", "mean_cov"]].sort_values("region_order")
    print(view.to_string(index=False))

if __name__ == "__main__":
    main()
