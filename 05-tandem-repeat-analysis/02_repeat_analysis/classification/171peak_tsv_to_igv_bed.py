#!/usr/bin/env python3
"""Convert 171peak_repeat_human.tsv to an IGV-viewable 9-column BED.

Follows the exact convention used in check_tandem_repeats_human.ipynb (cell 52/65):
  chrom       = chromosome
  start       = start - 1          (0-based)
  end         = end                (1-based, TRF end kept as-is)
  name        = TRF_<period_size>bp
  score       = alignment_score / max(alignment_score) * 1000  (int, 0-1000)
  strand      = .
  thickStart  = start - 1
  thickEnd    = end
  itemRgb     = color
"""
import csv
import sys

csv.field_size_limit(sys.maxsize)

TSV = "171peak_repeat_human.tsv"
BED = "171peak_repeat_human_igv.bed"


def main():
    rows = []
    with open(TSV) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)  # skip header
        for line in reader:
            rows.append(line)

    # Columns (0-indexed): 1=sequence 2=start 3=end 4=period_size 9=alignment_score
    #                       20=chromosome 24=color
    max_score = max(float(r[9]) for r in rows)

    out_lines = []
    for r in rows:
        chrom = r[20]
        start = int(r[2]) - 1  # 0-based
        end = int(r[3])
        period = int(float(r[4]))
        score = int(float(r[9]) / max_score * 1000)
        color = r[24]
        out_lines.append(
            f"{chrom}\t{start}\t{end}\tTRF_{period}bp\t{score}\t.\t{start}\t{end}\t{color}"
        )

    with open(BED, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")

    print(f"Wrote {len(out_lines)} intervals to {BED}")
    print(f"max alignment_score: {max_score:.1f}")


if __name__ == "__main__":
    main()
