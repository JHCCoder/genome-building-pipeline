#!/usr/bin/env python3
"""
Build pyCircularize_features_noMethylation.ipynb from the original
pyCircularize_features.ipynb: same circos figure minus the methylation track.

Data-prep comes verbatim from the original notebook's cells (fai, GC, GFF
features, A/B compartments, junctions). The plot cell is original cell 71 with:
  - methylation track block removed
  - methylation yticks block removed
  - A/B track re-spaced (55,45) -> (70,60)   (outer edge now flush with tRNA inner edge)
  - GC   track re-spaced (42,32) -> (57,47)  (shifted outward to preserve spacing)
  - latent bug fixed: max_dev_ab -> max_dev in ab_track.yticks
  - track_labels updated (drop CG methylation, restore 1)..8) numbering,
    move A/B & GC labels)

Writes a new notebook; does NOT touch the original.
"""
import json
import re

SRC = "pyCircularize_features.ipynb"
DST = "pyCircularize_features_noMethylation.ipynb"

nb = json.load(open(SRC))
cells = nb["cells"]
src = {i: "".join(cells[i]["source"]) for i in range(len(cells))}

# ---------------------------------------------------------------- data-prep
imports = src[1].replace("import pyBigWig\n", "")  # no longer needed

# fai: full read, all 30 chromosomes (cell 71 derives karyotype itself)
m = re.search(r'fai = pd\.read_csv\(\n(?:.*\n)*?\)\n', src[56])
assert m, "fai block not found"
fai_block = m.group(0)

# GC content
m = re.search(r'df_gc = pd\.read_csv\(\n(?:.*\n)*?\)\n', src[56])
assert m, "df_gc block not found"
gc_block = m.group(0)

# gene features: the GFF bound to seqid2genefeatures in cell 56 (peaks2utr)
lines = src[56].splitlines()
gff_line = next(l for l in lines if l.strip().startswith("gff_feat = Gff("))
seqid_line = next(l for l in lines if l.strip().startswith("seqid2genefeatures ="))
assert "peaks2utr" in gff_line, f"unexpected gff: {gff_line}"
gff_block = gff_line + "\n\n" + seqid_line + "\n"

# A/B compartments
m = re.search(r'df_ab = pd\.read_csv\(\n(?:.*\n)*?\)\n', src[56])
assert m, "df_ab block not found"
ab_block = m.group(0)

# junction bed
m = re.search(r'junction_bed = pd\.read_csv\(\n(?:.*\n)*?\)\n', src[56])
assert m, "junction_bed block not found"
junc_block = m.group(0)

# ---------------------------------------------------------------- plot cell
plot = src[71]

# 1) remove the whole methylation track block
m = re.search(
    r'    ####### Methylation track\n(?:.*\n)*?'
    r'    track\.fill_between\(\n(?:.*\n)*?    \)\n',
    plot, re.M)
assert m, "methylation block not found"
plot = plot.replace(m.group(0), "")

# 2) remove the methylation yticks block (references removed `track`)
m = re.search(
    r'        track\.yticks\(\n(?:.*\n)*?        \)\n',
    plot, re.M)
assert m, "methylation yticks block not found"
plot = plot.replace(m.group(0), "")

# 3) re-space A/B and GC (A/B outer edge = 70, flush with tRNA inner edge)
plot = plot.replace("ab_track = sector.add_track((55, 45))",
                    "ab_track = sector.add_track((70, 60))")
plot = plot.replace("gc_track = sector.add_track((42, 32))",
                    "gc_track = sector.add_track((57, 47))")

# 4) fix latent bug in ab_track.yticks labels
plot = plot.replace('labels=[f"{-max_dev_ab:.2f}","0", f"{max_dev_ab:.2f}"]',
                    'labels=[f"{-max_dev:.2f}","0", f"{max_dev:.2f}"]')

# 5) track_labels: drop methylation, restore numbering, move A/B & GC
plot = plot.replace('    (58.5, "CG methylation\\nfraction", "black"),\n', "")
plot = plot.replace('    (95, "Chromosomes", "black"),',
                    '    (95, "1) Chromosomes", "black"),')
plot = plot.replace('    (89.5, "Forward CDS", "tomato"),',
                    '    (89.5, "2) Forward CDS", "tomato"),')
plot = plot.replace('    (84.5, "Reverse CDS", "dodgerblue"),',
                    '    (84.5, "3) Reverse CDS", "dodgerblue"),')
plot = plot.replace('    (79.5, "snRNA", "limegreen"),',
                    '    (79.5, "4) snRNA", "limegreen"),')
plot = plot.replace('    (74.5, "lncRNA", "gold"),',
                    '    (74.5, "5) lncRNA", "gold"),')
plot = plot.replace('    (69.5, "tRNA", "mediumorchid"),',
                    '    (69.5, "6) tRNA", "mediumorchid"),')
plot = plot.replace('    (45.5, "A/B\\ncompartment", "black"),',
                    '    (60.5, "7) A/B\\ncompartment", "black"),')
plot = plot.replace('    (32.5, "GC\\ncontent", "black"),',
                    '    (47.5, "8) GC\\ncontent", "black"),')

# sanity: no residual methylation / bigWig references
for bad in ["bw_file", "bw.stats", "CG methylation", "max_dev_ab",
            "pyBigWig", "(68, 58)"]:
    assert bad not in plot, f"residual reference: {bad}"

# ---------------------------------------------------------------- save cells
save_png = '''fig.savefig(
    "degus_genome_circos_genetic_overview_noMethylation.png",
    dpi=600,
    bbox_inches="tight",
    transparent=False,
    facecolor="white"
)
'''
save_svg = '''fig.savefig(
    "degus_genome_circos_genetic_overview_noMethylation.svg",
    dpi=600,
    bbox_inches="tight",
    transparent=False,
    facecolor="white"
)
'''
save_pdf = '''fig.savefig(
    "degus_genome_circos_genetic_overview_noMethylation.pdf",
    dpi=600,
    bbox_inches="tight",
    transparent=False,
    facecolor="white"
)
'''

# ---------------------------------------------------------------- assemble
title_md = (
    "# Degus genome circos overview — no methylation track\n\n"
    "Same as `degus_genome_circos_genetic_epigenetic_overview` but the CpG "
    "methylation track is removed and the A/B compartment / GC content tracks "
    "are re-spaced outward to fill the freed band."
)


def code_cell(src_text):
    return {"cell_type": "code", "execution_count": None,
            "metadata": {}, "outputs": [],
            "source": src_text.splitlines(keepends=True)}


new_cells = [
    {"cell_type": "markdown", "metadata": {},
     "source": title_md.splitlines(keepends=True)},
    code_cell(imports),
    code_cell(fai_block),
    code_cell(gc_block),
    code_cell(gff_block),
    code_cell(ab_block),
    code_cell(junc_block),
    code_cell(plot),
    code_cell(save_png),
    code_cell(save_svg),
    code_cell(save_pdf),
]

out = {
    "cells": new_cells,
    "metadata": {
        "kernelspec": {
            "display_name": "python-visualizations",
            "language": "python",
            "name": "python-visualizations",
        },
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

with open(DST, "w") as f:
    json.dump(out, f, indent=1)

print(f"Wrote {DST}: {len(new_cells)} cells")
print("Plot cell edits OK; no methylation references remain.")
