#!/usr/bin/env python
"""
5a_matched_null.py -- generate matched-shuffle null placements for each
repeat array, matched on chromosome (exact), array length (exact), and window
covariates (mappability, GC content, repeat density) within tolerance.

A null that places homogeneous 349-bp arrays into more uniquely-mappable
non-repeat sequence would exaggerate apparent depletion; matching the shuffled
background on mappability / GC / repeat density prevents that.

Method: for each array, the "target" covariate vector is measured at the 100 kb
window containing the array's midpoint. All windows on the same chromosome
whose (mapp, GC, repd) all lie within tolerance of the target are the "allowed"
set. Null placements are then sampled uniformly from positions whose containing
window is in the allowed set (exact array length preserved). Arrays longer than
one window are handled naturally (their placement spans several windows, and the
containing-window check still applies at the placement start).

Output CSV columns:
  chrom family array_start array_end null_start null_end mapp gc repd matched
(matched=1 if >=1 allowed window existed; null placements only for matched)
"""
import argparse, random, csv, sys

def load_arrays(path, family):
    out = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            out.append((p[0], s, e, family))
    return out

def load_covar(path):
    # returns dict (chrom, start) -> value, and per-chrom sorted window starts
    d = {}
    starts = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                chrom, s, val = p[0], int(p[1]), float(p[3])
            except ValueError:
                continue
            d[(chrom, s)] = val
            starts.setdefault(chrom, []).append(s)
    for chrom in starts:
        starts[chrom].sort()
    return d, starts

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arrays-bin6", required=True)
    ap.add_argument("--arrays-bin4", required=True)
    ap.add_argument("--arrays-bin8", required=True)
    ap.add_argument("--mapp", required=True)
    ap.add_argument("--gc", required=True)
    ap.add_argument("--repd", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--n-shuf", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260807)
    ap.add_argument("--tol-mapp", type=float, default=0.10)
    ap.add_argument("--tol-gc", type=float, default=0.02)
    ap.add_argument("--tol-repd", type=float, default=0.10)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    W = 100000

    arrays = (load_arrays(args.arrays_bin6, "bin6_349") +
              load_arrays(args.arrays_bin4, "bin4_195") +
              load_arrays(args.arrays_bin8, "bin8_389"))

    mapp, mapp_starts = load_covar(args.mapp)
    gc,   gc_starts   = load_covar(args.gc)
    repd, repd_starts = load_covar(args.repd)
    all_starts = {c: sorted(set(mapp_starts[c]) & set(gc_starts[c]) & set(repd_starts[c]))
                  for c in mapp_starts}

    chrom_sizes = {}
    with open(args.chrom_sizes) as f:
        for line in f:
            p = line.split()
            if p[0].startswith("chr"):
                chrom_sizes[p[0]] = int(p[1])

    def covar_at(chrom, pos):
        s = (pos // W) * W
        return (mapp.get((chrom, s)), gc.get((chrom, s)), repd.get((chrom, s)))

    n_total = 0
    n_matched = 0
    n_null = 0
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chrom", "family", "array_start", "array_end",
                    "null_start", "null_end", "mapp", "gc", "repd", "matched"])
        for chrom, s, e, fam in arrays:
            if chrom not in all_starts or len(all_starts[chrom]) < 2:
                continue
            length = e - s
            clen = chrom_sizes[chrom]
            # target covariates at the array's own midpoint window
            tm, tg, tr = covar_at(chrom, (s + e) // 2)
            if tm is None or tg is None or tr is None:
                continue
            n_total += 1

            # allowed windows: same chromosome, covariate within tolerance
            allowed = []
            for ws in all_starts[chrom]:
                m_ = mapp.get((chrom, ws)); g_ = gc.get((chrom, ws)); r_ = repd.get((chrom, ws))
                if m_ is None or g_ is None or r_ is None:
                    continue
                if (abs(m_ - tm) <= args.tol_mapp and
                    abs(g_ - tg) <= args.tol_gc and
                    abs(r_ - tr) <= args.tol_repd):
                    allowed.append(ws)
            if not allowed:
                continue
            n_matched += 1

            # union of allowed window intervals -> allowed start range [lo,hi)
            intervals = []
            for ws in allowed:
                lo = ws
                hi = min(ws + W, clen)
                if lo < hi:
                    intervals.append((lo, hi))
            # merge intervals
            intervals.sort()
            merged = []
            for lo, hi in intervals:
                if merged and lo <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
                else:
                    merged.append((lo, hi))
            # feasible start range within [0, clen-length]
            feas = []
            for lo, hi in merged:
                s_lo = max(lo, 0)
                s_hi = min(hi, clen) - length
                if s_hi >= s_lo:
                    feas.append((s_lo, s_hi))
            if not feas:
                continue
            # total allowed length for uniform sampling
            total_len = sum(hi - lo + 1 for lo, hi in feas)

            for _ in range(args.n_shuf):
                # weighted pick among intervals
                r = rng.randint(0, total_len - 1)
                for lo, hi in feas:
                    span = hi - lo + 1
                    if r < span:
                        ns = lo + r
                        break
                    r -= span
                ne = ns + length
                w.writerow([chrom, fam, s, e, ns, ne,
                            f"{tm:.4f}", f"{tg:.4f}", f"{tr:.4f}", 1])
                n_null += 1

    sys.stderr.write(f"arrays processed: {n_total}, matched: {n_matched}, null placements: {n_null}\n")

if __name__ == "__main__":
    main()
