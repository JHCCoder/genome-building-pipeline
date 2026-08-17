#!/usr/bin/env python
# ============================================================================
# 7_methylation_cpg_metrics.py -- CpG callability + mappability metrics around
# the peak-based CENP-A cores, in the same spatial intervals as the occupancy
# analysis. Answers: is the sparse core signal due to (a) few CpG positions in
# the reference, (b) CpGs present but not callable (no reads / unmappable)?
#
# Per interval (chrom x region) metrics:
#   ref_CG            # CG dinucleotides in the reference assembly sequence
#   cg_density_per_kb = 1000 * ref_CG / interval_bp        (potential CpG sites)
#   n_cov0 / n_cov5   # CpGs with coverage >0 / >=MIN_COV   (callable)
#   callable_frac     n_cov5 / ref_CG
#   mean_frac, wmean_frac   methylation fraction at callable CpGs
#   n_meth5, meth_density_per_kb   callable & frac>=0.5
#   mean_cov          read depth at covered CpGs
#   mappability       fraction of interval bp covered by pooled CUT&Tag
#                     control (k=1) fragments, aggregated from 100 kb windows
#                     (data/mappability/mappability_win100kb.bed) -- the same
#                     mappability covariate used by the matched-null analysis
#
# Usage: python 7_methylation_cpg_metrics.py <work_dir>
# ============================================================================
import os, sys, subprocess
import numpy as np
import pandas as pd
import pyBigWig

WORK = sys.argv[1] if len(sys.argv) > 1 else \
    "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-repeat-chromosome"
DOM    = os.path.join(WORK, "data", "domains", "cenpa_domains_peakbased.csv")
MAP    = os.path.join(WORK, "data", "mappability", "mappability_win100kb.bed")
BW_DIR = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/circos-plot/feature-overview"
FRAC   = os.path.join(BW_DIR, "degu_6834_PFC_1.CGN-both.frac.bw")
COV    = os.path.join(BW_DIR, "degu_6834_PFC_1.CGN-both.cov.bw")
FASTA  = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm-041425-assembly-mitoFiltered-scaffolded-curated-masked-chrNameAssigned-contamFiltered/assembly_final.sorted.headerRenamed.chrAssigned.contamFiltered.fasta"
SAM    = "/tscc/projects/ps-renlab2/jhc103/toolshed/cactus-bin-v2.9.2/bin/samtools"
OUT    = os.path.join(WORK, "results", "methylation_cpg_metrics.csv")

MIN_COV = 5
METH_CUT = 0.5

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

def chrom_seq(chrom):
    """Fetch full chromosome sequence as bytes (uppercase+lowercase preserved)."""
    r = subprocess.run([SAM, "faidx", FASTA, chrom], capture_output=True)
    if r.returncode != 0:
        return None
    return b"".join(r.stdout.splitlines()[1:])

def count_cg(seq, s, e):
    """Case-insensitive CG dinucleotide count over seq[s:e]."""
    if e > len(seq): e = len(seq)
    sub = seq[s:e]
    return sub.count(b"CG") + sub.count(b"cg")

def main():
    dom = pd.read_csv(DOM)
    fr = pyBigWig.open(FRAC); cv = pyBigWig.open(COV)

    # mappability windows -> (chrom, start, end, frac)
    mw = pd.read_csv(MAP, sep="\t", header=None, names=["chrom", "start", "end", "frac"])
    mw_g = mw.groupby("chrom")

    cache = {}
    rows = []
    for _, d in dom.iterrows():
        chrom = str(d["chrom"]); cs = int(d["core_start"]); ce = int(d["core_end"])
        clen = int(d["chrom_len"])
        seq = cache.get(chrom)
        if seq is None:
            seq = chrom_seq(chrom); cache[chrom] = seq
        if seq is None:
            print(f"  !! samtools failed for {chrom} -- skipping CG counts"); seq = b""

        for (side, region, s, e) in build_intervals(cs, ce, clen):
            key = side if side in ("core", "background") else f"{side}_{region}"
            ibp = e - s

            # reference CpG sites
            ref_cg = count_cg(seq, s, e) if seq else np.nan

            # callable CpGs from the WGBS bigwigs
            fi = fr.intervals(chrom, s, e) or []
            ci = cv.intervals(chrom, s, e) or []
            covd = {st: v for (st, _, v) in ci}
            fracs, covs = [], []
            for (st, en, fv) in fi:
                c = covd.get(st, np.nan)
                if not np.isnan(c) and c > 0:
                    fracs.append(fv); covs.append(c)
            n_cov0 = len(fracs)
            sel = [i for i, c in enumerate(covs) if c >= MIN_COV]
            n_cov5 = len(sel)
            f5 = [fracs[i] for i in sel]; c5 = [covs[i] for i in sel]
            n_meth5 = sum(1 for f in f5 if f >= METH_CUT)
            mean_frac = float(np.mean(f5)) if f5 else np.nan
            mean_cov  = float(np.mean(c5)) if c5 else np.nan

            # mappability (100 kb window aggregate, weighted by overlap bp)
            w = mw_g.get_group(chrom) if chrom in mw_g.groups else pd.DataFrame()
            if len(w):
                ov = np.minimum(w["end"], e) - np.maximum(w["start"], s)
                ov = np.maximum(ov, 0)
                mapp = float((ov * w["frac"]).sum()) / ibp
            else:
                mapp = np.nan

            rows.append({
                "chromosome": chrom, "core_start": cs, "core_end": ce,
                "side": side, "region": region, "region_order": REGION_ORDER[key],
                "interval_start": s, "interval_end": e, "interval_bp": ibp,
                "ref_CG": int(ref_cg), "cg_density_per_kb": 1000 * ref_cg / ibp,
                "n_cov0": n_cov0, "n_cov5": n_cov5,
                "callable_frac": n_cov5 / ref_cg if ref_cg else np.nan,
                "n_meth5": n_meth5,
                "meth_density_per_kb": 1000 * n_meth5 / ibp,
                "mean_frac": mean_frac, "mean_cov": mean_cov,
                "mappability": mapp,
            })
    fr.close(); cv.close()

    met = pd.DataFrame(rows)
    met.to_csv(OUT, index=False)

    # genome-wide pooled rows (weighted by n_cov5 across chromosomes)
    gw = []
    for (side, region, ro), g in met.groupby(["side", "region", "region_order"]):
        g = g[g["n_cov5"] > 0]
        ibp = int(g["interval_bp"].sum())
        ref = int(g["ref_CG"].sum())
        n0, n5 = int(g["n_cov0"].sum()), int(g["n_cov5"].sum())
        nm5 = int(g["n_meth5"].sum())
        row = {"chromosome": "genome", "side": side, "region": region,
               "region_order": ro, "interval_bp": ibp,
               "ref_CG": ref, "cg_density_per_kb": 1000 * ref / ibp,
               "n_cov0": n0, "n_cov5": n5,
               "callable_frac": n5 / ref if ref else np.nan, "n_meth5": nm5,
               "meth_density_per_kb": 1000 * nm5 / ibp,
               "mean_frac": float(np.average(g["mean_frac"], weights=g["n_cov5"])),
               "mean_cov": float(np.average(g["mean_cov"], weights=g["n_cov5"])),
               "mappability": float(np.average(g["mappability"], weights=g["interval_bp"]))}
        gw.append(row)
    met = pd.concat([met, pd.DataFrame(gw)], ignore_index=True)
    met.to_csv(OUT, index=False)

    print(f"wrote {OUT}")
    view = met[met["chromosome"] == "genome"][
        ["region_order", "side", "region", "interval_bp", "ref_CG",
         "cg_density_per_kb", "n_cov0", "n_cov5", "callable_frac",
         "n_meth5", "meth_density_per_kb", "mean_frac", "mean_cov", "mappability"]
    ].sort_values("region_order")
    with pd.option_context("display.width", 250):
        print(view.to_string(index=False))

if __name__ == "__main__":
    main()
