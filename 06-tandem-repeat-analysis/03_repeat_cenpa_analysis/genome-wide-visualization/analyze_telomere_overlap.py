#!/usr/bin/env python3
"""
Compare telomeric (TTAGGG/CCCTAA) arrays against:
  - 349 bp satellite arrays (bin6)
  - centroAnno-predicted centromeric repeat regions
  - CENP-A positive intervals

For each 349 bp array we report:
  - length
  - whether it overlaps a telomere array (and how much)
  - distance to nearest telomere array
  - distance to nearest chromosome end
  - whether it overlaps centroAnno / CENP-A

Output: summary table + per-array table.
"""
import argparse, bisect, sys


def load_bed(path):
    out = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue  # header row
            out.setdefault(p[0], []).append((s, e))
    for c in out:
        out[c].sort()
    return out


def merge_intervals(ivs):
    ivs = sorted(ivs)
    out = []
    for s, e in ivs:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def overlap_length(a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    return max(0, hi - lo)


def nearest_distance(ivs, s, e):
    """min distance from [s,e) to any interval in a sorted list of (start,end).
    Uses interval END for intervals to the left, START for intervals to the right.
    Overlap/adjacency -> 0."""
    starts = [iv[0] for iv in ivs]
    i = bisect.bisect_right(starts, s) - 1
    best = None
    if i >= 0:
        # interval to the left: distance = s - end of that interval
        d = s - ivs[i][1]
        best = max(0, d)
    if i + 1 < len(ivs):
        d = ivs[i + 1][0] - e
        best = d if best is None else min(best, max(0, d))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--telomere", required=True)
    ap.add_argument("--arrays", required=True)
    ap.add_argument("--centro", required=True)
    ap.add_argument("--cenpa", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--min-telomere-bp", type=int, default=50,
                    help="ignore telomere arrays shorter than this")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    # Chromosome sizes
    chrom_sizes = {}
    with open(args.fai) as f:
        for line in f:
            p = line.split("\t")
            chrom_sizes[p[0]] = int(p[1])

    # Telomere arrays -> per chromosome list of (s,e), then merge same-strand-agnostic
    tel = load_bed(args.telomere)
    tel = {c: [iv for iv in ivs if iv[1] - iv[0] >= args.min_telomere_bp]
           for c, ivs in tel.items()}
    tel = {c: merge_intervals(ivs) for c, ivs in tel.items() if ivs}

    centro = load_bed(args.centro)
    cenpa = load_bed(args.cenpa)
    arrays = load_bed(args.arrays)

    n_tel = sum(len(v) for v in tel.values())
    tel_bp = sum(e - s for v in tel.values() for s, e in v)
    print(f"Telomere arrays (>= {args.min_telomere_bp} bp): {n_tel}; total {tel_bp/1e3:.1f} kb",
          file=sys.stderr)

    # Per-array table
    rows = []
    for c in sorted(arrays):
        clen = chrom_sizes.get(c, 0)
        tel_ivs = tel.get(c, [])
        for s, e in arrays[c]:
            ln = e - s
            # overlap with telomere
            tel_overlap = 0
            for ts, te in tel_ivs:
                tel_overlap += overlap_length((s, e), (ts, te))
            dist_tel = nearest_distance(tel_ivs, s, e) if tel_ivs else None
            dist_start = s
            dist_end = (clen - e) if clen else None
            at_end = (dist_start <= 200000) or (dist_end is not None and dist_end <= 200000)
            ov_centro = sum(overlap_length((s, e), (cs, ce)) for cs, ce in centro.get(c, []))
            ov_cenpa = sum(overlap_length((s, e), (cs, ce)) for cs, ce in cenpa.get(c, []))
            rows.append((c, s, e, ln, tel_overlap, dist_tel, dist_start, dist_end, at_end, ov_centro, ov_cenpa))

    with open(args.out_prefix + "_per_array.tsv", "w") as f:
        f.write("chrom\tstart\tend\tlen\tov_telomere_bp\tdist_telomere_bp\tdist_start\tdist_end\tat_chr_end\toverlap_centro_bp\toverlap_cenpa_bp\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    # Summary
    n = len(rows)
    ov_tel = [r for r in rows if r[4] > 0]
    # arrays that are "at chromosome end" and check proximity to telomere
    end_arrays = [r for r in rows if r[8]]
    near_tel = [r for r in rows if r[5] is not None and r[5] <= 50000]
    print("\n=== Summary ===")
    print(f"349bp arrays: {n}")
    print(f"  overlap TTAGGG array (bp>0): {len(ov_tel)} ({100*len(ov_tel)/n:.1f}%)")
    print(f"  within 50 kb of a TTAGGG array: {len(near_tel)} ({100*len(near_tel)/n:.1f}%)")
    print(f"  at chromosome end (<=200 kb): {len(end_arrays)} ({100*len(end_arrays)/n:.1f}%)")
    if end_arrays:
        end_near_tel = [r for r in end_arrays if r[5] is not None and r[5] <= 50000]
        end_far = [r for r in end_arrays if r[5] is None or r[5] > 50000]
        print(f"    of end arrays: within 50 kb of TTAGGG {len(end_near_tel)}; far/absent {len(end_far)}")
    # overlap centro/cenpa
    ov_c = [r for r in rows if r[9] > 0]
    ov_p = [r for r in rows if r[10] > 0]
    print(f"  overlap centroAnno: {len(ov_c)} ({100*len(ov_c)/n:.1f}%)")
    print(f"  overlap CENP-A positive: {len(ov_p)} ({100*len(ov_p)/n:.1f}%)")
    # lengths of end vs internal
    end_bp = sum(r[3] for r in end_arrays)
    int_bp = sum(r[3] for r in rows if not r[8])
    print(f"  bp at chromosome ends: {end_bp/1e6:.2f} Mb ({100*end_bp/(end_bp+int_bp):.1f}%)")


if __name__ == "__main__":
    main()
