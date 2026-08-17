#!/usr/bin/env python3
"""Convert classifier_3367_caseinsensitive.tsv -> IGV-viewable 9-column BEDs.

Reads the SW-classification of all 3,367 period-110-135 TRF hits and writes:

  mouse_period120_classified_igv.bed         -- all hits, colored by class
  mouse_minor_satellite_confirmed_igv.bed    -- ONLY confirmed minor-satellite
                                                hits (the genuine ones)

Colors:
  minor      #e41a1c  (red)     -- confirmed minor satellite
  major      #ff7f00  (orange)  -- major satellite (pericentromere)
  ambiguous  #984ea3  (purple)  -- matches both, cannot separate
  other      #d3d3d3  (light grey) -- TRF period false positive (not satellite)
"""
import csv

CWD = "/tscc/nfs/home/jhc103/ps-renlab2-link/degu-genome-assembly-proj/code/command-line-script/genome-annotation/trf-tandem-repeat"
TSV = CWD + "/classifier_3367_caseinsensitive.tsv"

COLOR = {
    "minor": "#e41a1c",
    "major": "#ff7f00",
    "ambiguous": "#984ea3",
    "other": "#d3d3d3",
    "NO_SEQ": "#999999",
}
SCORE = {
    "minor": 1000,
    "major": 800,
    "ambiguous": 600,
    "other": 200,
    "NO_SEQ": 100,
}


def main():
    rows = []
    with open(TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for x in r:
            rows.append(x)
    print(f"read {len(rows)} rows")

    all_lines = []
    minor_lines = []
    for x in rows:
        cat = x["class"]
        chrom, s, e = x["chrom"], int(x["start"]), int(x["end"])
        score = SCORE[cat]
        color = COLOR[cat]
        smi = x["smi"]
        name = f"{cat}_smi{smi}"
        all_lines.append(f"{chrom}\t{s}\t{e}\t{name}\t{score}\t.\t{s}\t{e}\t{color}")
        if cat == "minor":
            minor_lines.append(f"{chrom}\t{s}\t{e}\t{name}\t1000\t.\t{s}\t{e}\t{color}")

    with open(CWD + "/mouse_period120_classified_igv.bed", "w") as fh:
        fh.write("\n".join(all_lines) + "\n")
    with open(CWD + "/mouse_minor_satellite_confirmed_igv.bed", "w") as fh:
        fh.write("\n".join(minor_lines) + "\n")

    from collections import Counter
    print("class counts:", dict(Counter(x["class"] for x in rows)))
    print(f"wrote mouse_period120_classified_igv.bed ({len(all_lines)} intervals)")
    print(f"wrote mouse_minor_satellite_confirmed_igv.bed ({len(minor_lines)} intervals)")


if __name__ == "__main__":
    main()
