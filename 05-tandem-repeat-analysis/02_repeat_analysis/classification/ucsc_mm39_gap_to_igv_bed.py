#!/usr/bin/env python3
"""Extract UCSC mm39 gap table -> IGV-viewable BED, plus cross-check vs minor satellite.

Official source (downloaded 2026-08-06):
  https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/gap.txt.gz
  Columns: bin, chrom, chromStart, chromEnd, ix, n, size, type, bridged
  (chromStart/chromEnd are already 0-based half-open.)

UCSC chr names are mapped to RefSeq accessions (NC_XXXXXXXX.N) using the
local GRCm39 FASTA headers, so the BEDs load against the NCBI FASTA in IGV.

Outputs (all 9-column IGV BED):
  mm39_UCSC_gap_centromeres_igv.bed        -- 20 centromere gaps, NC-named
  mm39_UCSC_gap_centromeres_chrName.bed    -- 20 centromere gaps, chr-named
  mm39_UCSC_gap_all_igv.bed                -- all 324 gaps, NC-named, colored by type
"""
import gzip
import re
import sys

cwd = ("/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/"
       "code/command-line-script/genome-annotation/trf-tandem-repeat")
GAP = "/tscc/lustre/ddn/scratch/jhc103/ucsc_mm39_gap/mm39_gap.txt.gz"
FASTA = ("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/"
         "data/GRCm39_genome/GCF_000001635.27_GRCm39_genomic.fna")
MINOR_BED = cwd + "/mouse_satellite_minor120_origChromName_igv.bed"

COLOR = {
    "centromere": "#e41a1c",
    "telomere":   "#377eb8",
    "short_arm":  "#4daf4a",
    "contig":     "#984ea3",
    "scaffold":   "#999999",
}


def load_chr_map(fasta):
    """UCSC chrN -> RefSeq accession, parsed from FASTA headers."""
    m = {}
    with open(fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            acc = line[1:].split()[0]
            mm = re.search(r"chromosome (\d+|X|Y),", line)
            if mm:
                m[f"chr{mm.group(1)}"] = acc
    return m


def load_gaps(path):
    """Return list of (chrom, start, end, size, type)."""
    gaps = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gaps.append((f[1], int(f[2]), int(f[3]), int(f[6]), f[7]))
    return gaps


def load_bed_intervals(path):
    """Return dict chrom -> list of (start, end)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return out


def bed_line(chrom, start, end, gtype):
    name = f"UCSC_gap_{gtype}"
    color = COLOR.get(gtype, "#000000")
    score = 1000 if gtype == "centromere" else 500
    return f"{chrom}\t{start}\t{end}\t{name}\t{score}\t.\t{start}\t{end}\t{color}"


def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def main():
    chr_map = load_chr_map(FASTA)
    gaps = load_gaps(GAP)

    # ---- write BEDs -----------------------------------------------------
    cent = [g for g in gaps if g[4] == "centromere"]
    print(f"gaps total: {len(gaps)}  |  centromere gaps: {len(cent)}")

    # NC-named centromere BED
    with open(cwd + "/mm39_UCSC_gap_centromeres_igv.bed", "w") as fh:
        for chrom, s, e, size, gtype in cent:
            fh.write(bed_line(chr_map.get(chrom, chrom), s, e, gtype) + "\n")
    print(f"wrote mm39_UCSC_gap_centromeres_igv.bed ({len(cent)} intervals, NC-named)")

    # chr-named centromere BED
    with open(cwd + "/mm39_UCSC_gap_centromeres_chrName.bed", "w") as fh:
        for chrom, s, e, size, gtype in cent:
            fh.write(bed_line(chrom, s, e, gtype) + "\n")
    print(f"wrote mm39_UCSC_gap_centromeres_chrName.bed ({len(cent)} intervals, chr-named)")

    # all gaps, NC-named
    with open(cwd + "/mm39_UCSC_gap_all_igv.bed", "w") as fh:
        for chrom, s, e, size, gtype in gaps:
            fh.write(bed_line(chr_map.get(chrom, chrom), s, e, gtype) + "\n")
    print(f"wrote mm39_UCSC_gap_all_igv.bed ({len(gaps)} intervals, NC-named)")

    # ---- cross-check vs minor satellite ----------------------------------
    minor = load_bed_intervals(MINOR_BED)
    chrom_of_gap = {}
    print("\nUCSC centromere gap vs minor-satellite TRF intervals (NC-named):")
    print(f"{'chrom':<14}{'NC accession':<16}{'gap (Mb)':<24}{'# minor ovlp':<14}{'max ovlp bp':<12}")
    for chrom, s, e, size, gtype in sorted(cent, key=lambda g: g[0]):
        nc = chr_map.get(chrom, chrom)
        ivs = minor.get(nc, [])
        cnt = 0
        max_ov = 0
        for (ms, me) in ivs:
            ov = overlap(ms, me, s, e)
            if ov > 0:
                cnt += 1
                max_ov = max(max_ov, ov)
        print(f"{chrom:<14}{nc:<16}{s/1e6:.2f}-{e/1e6:.2f}{cnt:<14}{max_ov:<12}")

    # chrY: no UCSC centromere gap, but RefSeq GFF annotates one
    y_ivs = minor.get("NC_000087.8", [])
    y_ov = sum(overlap(ms, me, 4072168, 4161965) for ms, me in y_ivs)
    y_cnt = sum(1 for ms, me in y_ivs if overlap(ms, me, 4072168, 4161965) > 0)
    print("\nchrY (NC_000087.8) RefSeq GFF centromere 4,072,168-4,161,965:")
    print(f"  UCSC gap table centromere gap: NONE")
    print(f"  minor-satellite intervals overlapping GFF centromere: {y_cnt}  (total overlap {y_ov} bp)")


if __name__ == "__main__":
    main()
