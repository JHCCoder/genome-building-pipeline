#!/usr/bin/env python3
"""prepare_panel_data.py — build the 4 cached data tables for the
"195 bp repeats are L1" 3-panel supplementary figure.

  panel_a_composition.tsv  — RepeatMasker class composition of bin4 (195 bp)
                             vs bin6 (349 bp) arrays as an EXCLUSIVE partition:
                             each base assigned to exactly one category, with
                             priority L1(fam-189 / fam-18 / other) > Unknown
                             > Other > Unannotated. Bars therefore sum to 100%.
  panel_b_consensus.tsv    — L1 consensus 5' 195 bp tandem-repeat blocks
                             (from TRF on the consensus, already verified).
  panel_c_locus.tsv        — RepeatMasker annotations across the example locus
                             chr21:83,255,000-83,340,000 (the longest 195 bp
                             array), colored by family/class.
  panel_c_arrays.tsv       — bin4 (195 bp) arrays in that window.

The plot-only R script reads these tables; re-run this only when data changes.
"""
import os, sys

WD = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/figure/cenpa-cuttag-enrichment/195bp-repeatmasker-overlap"
SC = "/tscc/lustre/ddn/scratch/jhc103"
FAM = os.path.join(SC, "rm_hifiasm041425.family.bed")
RES = os.path.join(WD, "results")
os.makedirs(RES, exist_ok=True)

def load_bed(path, ncols):
    out = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            chrom, s, e = p[0], int(p[1]), int(p[2])
            out.append((chrom, s, e) + tuple(p[3:]))
    return out

def union_coverage(cat, arrays):
    """bp of category elements (chrom,start,end) covered by array intervals,
    counted once. arrays: list of (chrom,start,end)."""
    # clip category elements to arrays, then merge per chromosome
    from collections import defaultdict
    clipped = defaultdict(list)
    # sort arrays by chrom,start
    arr = defaultdict(list)
    for a in arrays:
        arr[a[0]].append((a[1], a[2]))
    for chrom, segs in arr.items():
        segs.sort()
    # for each category element, find overlapping array segments (merge first per chrom)
    cat = defaultdict(list)
    for c in cat:
        cat[c[0]].append((c[1], c[2]))
    for chrom, segs in cat.items():
        segs.sort()
        if chrom not in arr:
            continue
        asegs = arr[chrom]
        i = 0
        for (s, e) in segs:
            # advance asegs pointer
            while i < len(asegs) and asegs[i][1] <= s:
                i += 1
            j = i
            while j < len(asegs) and asegs[j][0] < e:
                a = max(s, asegs[j][0]); b = min(e, asegs[j][1])
                if b > a:
                    clipped[chrom].append((a, b))
                j += 1
    # merge clipped per chrom
    total = 0
    for chrom, segs in clipped.items():
        segs.sort()
        for (s, e) in segs:
            total += e - s
    return total

# ── 1. Load category element files ──────────────────────────────────────────
cats = {}
for name in ("l1_fam189", "l1_fam18", "l1_other", "unknown", "other"):
    p = os.path.join(SC, "rmcats", name + ".bed")
    cats[name] = load_bed(p, 3)

# ── 2. Load arrays ──────────────────────────────────────────────────────────
arrays = {}
for lbl, fn in (("bin4_195bp", "bin4_arrays_seq.bed"), ("bin6_349bp", "bin6_arrays_seq.bed")):
    arrays[lbl] = load_bed(os.path.join(WD, fn), 5)

# ── 3. Exclusive partition (priority L1 > Unknown > Other > Unannotated) ────
# Work per chromosome with a sorted sweep. L1 union first.
def exclusive_partition(arr_intervals, cats):
    # returns dict category->bp
    from collections import defaultdict
    # build per-chrom sorted array segments
    arr_by = defaultdict(list)
    for a in arr_intervals:
        arr_by[a[0]].append((a[1], a[2]))
    for chrom in arr_by:
        arr_by[chrom].sort()
    result = defaultdict(int)
    # 1. L1 coverage (union of fam189+fam18+other)
    l1_union = defaultdict(list)
    for name in ("l1_fam189", "l1_fam18", "l1_other"):
        for c in cats[name]:
            l1_union[c[0]].append((c[1], c[2]))
    for chrom in l1_union:
        l1_union[chrom].sort()
    # clip L1 to arrays
    l1_clip = defaultdict(list)
    for chrom, lsegs in l1_union.items():
        if chrom not in arr_by:
            continue
        ase = arr_by[chrom]; i = 0
        for (s, e) in lsegs:
            while i < len(ase) and ase[i][1] <= s:
                i += 1
            j = i
            while j < len(ase) and ase[j][0] < e:
                a = max(s, ase[j][0]); b = min(e, ase[j][1])
                if b > a:
                    l1_clip[chrom].append((a, b))
                j += 1
    # merge l1_clip, and count L1 per family + total
    l1_fam = {"l1_fam189": defaultdict(list), "l1_fam18": defaultdict(list), "l1_other": defaultdict(list)}
    for name in ("l1_fam189", "l1_fam18", "l1_other"):
        for c in cats[name]:
            l1_fam[name][c[0]].append((c[1], c[2]))
    for name in l1_fam:
        for chrom in l1_fam[name]:
            l1_fam[name][chrom].sort()
    # per-family coverage within arrays
    for name in ("l1_fam189", "l1_fam18", "l1_other"):
        bp = 0
        for chrom, segs in l1_fam[name].items():
            if chrom not in arr_by:
                continue
            ase = arr_by[chrom]; i = 0
            for (s, e) in segs:
                while i < len(ase) and ase[i][1] <= s:
                    i += 1
                j = i
                while j < len(ase) and ase[j][0] < e:
                    a = max(s, ase[j][0]); b = min(e, ase[j][1])
                    if b > a:
                        bp += b - a
                    j += 1
        result[name] = bp
    # 2. non-L1 remainder of arrays
    def subtract(arr_segs, sub_segs):
        out = []
        i = 0
        for (s, e) in arr_segs:
            cur = s
            while i < len(sub_segs) and sub_segs[i][1] <= cur:
                i += 1
            j = i
            while j < len(sub_segs) and sub_segs[j][0] < e:
                if sub_segs[j][0] > cur:
                    out.append((cur, sub_segs[j][0]))
                cur = max(cur, sub_segs[j][1])
                j += 1
            if cur < e:
                out.append((cur, e))
        return out
    nonl1 = {}
    for chrom, segs in arr_by.items():
        srt = sorted(l1_clip[chrom])
        nonl1[chrom] = subtract(segs, srt)
    # 3. Unknown and Other within the non-L1 remainder
    def cov_within(cat_segs_by_chrom, region_by_chrom):
        bp = 0
        for chrom, segs in cat_segs_by_chrom.items():
            reg = sorted(region_by_chrom.get(chrom, []))
            if not reg:
                continue
            reg.sort()
            i = 0
            for (s, e) in segs:
                while i < len(reg) and reg[i][1] <= s:
                    i += 1
                j = i
                while j < len(reg) and reg[j][0] < e:
                    a = max(s, reg[j][0]); b = min(e, reg[j][1])
                    if b > a:
                        bp += b - a
                    j += 1
        return bp
    unk_segs = defaultdict(list)
    for c in cats["unknown"]:
        unk_segs[c[0]].append((c[1], c[2]))
    for chrom in unk_segs:
        unk_segs[chrom].sort()
    oth_segs = defaultdict(list)
    for c in cats["other"]:
        oth_segs[c[0]].append((c[1], c[2]))
    for chrom in oth_segs:
        oth_segs[chrom].sort()
    result["unknown"] = cov_within(unk_segs, nonl1)
    result["other"] = cov_within(oth_segs, nonl1)
    # 4. unannotated
    total = sum(e - s for segs in arr_by.values() for s, e in segs)
    result["unannotated"] = max(0, total - sum(result.values()))
    return result, total

order = ["l1_fam189", "l1_fam18", "l1_other", "unknown", "other", "unannotated"]
with open(os.path.join(RES, "panel_a_composition.tsv"), "w") as f:
    f.write("array_group\tcategory\tbp\tpct\n")
    for lbl in ("bin4_195bp", "bin6_349bp"):
        part, total = exclusive_partition(arrays[lbl], cats)
        for cat in order:
            bp = part[cat]
            f.write(f"{lbl}\t{cat}\t{bp}\t{100*bp/total:.4f}\n")
print("panel_a written:")
print(open(os.path.join(RES, "panel_a_composition.tsv")).read())

# ── 4. Panel (b): L1 consensus 5' repeat blocks (from verified TRF) ─────────
with open(os.path.join(RES, "panel_b_consensus.tsv"), "w") as f:
    f.write("family\tconsensus_len\trepeat_start\trepeat_end\tunit_len\tblocks_csv\n")
    f.write("rnd-1_family-189\t7758\t1\t804\t195\t1-195|196-390|391-585|586-780|781-804\n")
    f.write("rnd-1_family-18\t4933\t1\t438\t195\t1-195|196-390|391-438\n")
print("panel_b written")

# ── 5. Panel (c): example locus chr21 (seq22):83,255,000-83,340,000 ────────
W0, W1 = 83255000, 83340000
def fam_color(fam, cls):
    if cls == "LINE/L1":
        if fam == "rnd-1_family-189": return "L1 fam-189"
        if fam == "rnd-1_family-18":  return "L1 fam-18"
        return "L1 other"
    if cls == "Unknown": return "Unknown"
    return "Other"
rows = []
with open(FAM) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6: continue
        chrom, s, e, fam, cls, strand = p[0], int(p[1]), int(p[2]), p[3], p[4], p[5]
        if chrom == "seq22" and s < W1 and e > W0:
            rows.append((s, e, fam_color(fam, cls), fam, strand))
rows.sort()
with open(os.path.join(RES, "panel_c_locus.tsv"), "w") as f:
    f.write("chr\tstart\tend\tcategory\tfamily\tstrand\n")
    for s, e, c, fam, strand in rows:
        f.write(f"chr21\t{s}\t{e}\t{c}\t{fam}\t{strand}\n")
print(f"panel_c_locus written: {len(rows)} RM annotations in window")
# bin4 arrays in the window
arrs = [a for a in arrays["bin4_195bp"] if a[0] == "seq22" and a[1] < W1 and a[2] > W0]
with open(os.path.join(RES, "panel_c_arrays.tsv"), "w") as f:
    f.write("chr\tstart\tend\tarray_id\n")
    for a in sorted(arrs, key=lambda x: x[1]):
        f.write(f"chr21\t{a[1]}\t{a[2]}\t{a[3]}\n")
print(f"panel_c_arrays written: {len(arrs)} bin4 arrays in window")
