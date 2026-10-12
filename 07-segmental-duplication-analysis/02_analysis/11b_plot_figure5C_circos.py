#!/usr/bin/env python3
"""
Figure 5C: genome-wide circos of SegDup hotspots and the interferon-response genes.

Same tracks and data as the second plotting cell of 11_circos_segdup_hotspot.ipynb
(pyCirclize), redrawn at printed size for Genome Research. From outside in:
  1. chromosomes 1-28, X, Y (ticks every 50 Mb, numbers outside);
  2. SegDup-hotspot density: % of each 5-Mb window (100-kb step) covered by the
     1,425 hotspot regions (segdup_density_rolling_..._hotspot80perc.tsv);
  3. the 1,271 genes in SegDup hotspots (grey ticks);
  4. the eight GO:0071346 genes (cellular response to type II interferon) and their
     paralog copies, coloured by family (members from merged_paralog_singlecopy_families.tsv);
  5. SegDup links (BISER) that overlap one of those genes, coloured by family, drawn as
     stroked curves between region midpoints (no filled ribbons; see note at circos.link_line).
Track descriptions live in the figure legend; the in-panel numbered headings were removed.
The family key below the circle uses italic gene symbols (mouse-style nomenclature).
Colours: Paul Tol "muted" (colour-blind safe); density in dark grey.

Run with the python-visualizations env (pyCirclize 1.10):
    python 11b_plot_figure5C_circos.py [--width-mm 90] [--height-mm 100]
"""

import argparse
import os
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqFeature import FeatureLocation, SeqFeature
from pycirclize import Circos


def _cfg(name):
    val = os.environ.get(name)
    if val:
        return val
    root = os.path.abspath(os.path.dirname(__file__))
    while not os.path.isfile(os.path.join(root, "config.sh")) and root != os.path.dirname(root):
        root = os.path.dirname(root)
    # Source config.sh in a subshell so nested $VAR references are expanded.
    val = subprocess.check_output(
        ["bash", "-c", f'source "$1" 2>/dev/null; printf "%s" "${{{name}}}"', "_",
         os.path.join(root, "config.sh")],
        text=True,
    )
    if not val:
        sys.exit(f"ERROR: {name} is not set in the environment or in config.sh")
    return val


FONT_SIZE = 8
GENE_COLORS = {  # gene family -> colour (Paul Tol muted)
    "Gbp1": "#332288", "Gbp5": "#117733", "Gbp6": "#AA4499", "Ifng": "#CC6677",
    "Casp1": "#DDCC77", "Ccl3": "#44AA99", "Lgals9": "#882255", "Syncrip": "#88CCEE",
}
DENSITY_COLOR = "#5A5A5A"
HOTSPOT_GENE_COLOR = "#A8A8A8"

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=90)
ap.add_argument("--height-mm", type=float, default=100)
opts = ap.parse_args()
PANEL_W, PANEL_H = opts.width_mm / 25.4, opts.height_mm / 25.4
here = os.path.dirname(os.path.abspath(__file__))
out_dir = _cfg("FIG5_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, "figure5C_circos")

# ---------------------------------------------------------------- data
fai = pd.read_csv(_cfg("FIG5_FAI"), sep="\t", header=None, usecols=[0, 1], names=["chr", "length"])
karyotype = fai[fai["chr"].str.contains("chr")].head(30)
chr_sizes = dict(zip(karyotype["chr"], karyotype["length"]))

density = pd.read_csv(os.path.join(here, "segdup_density_rolling_5000000_100000_bpLevel_hotspot80perc.tsv"),
                      sep="\t", index_col=False)
hot_genes = pd.read_csv(os.path.join(here, "segdup_genes_hotspots_80percLength.tsv"), sep="\t", index_col=False)
print(f"hotspot genes drawn: {hot_genes['gene'].nunique():,}")

bedpe = pd.read_csv(_cfg("FIG5_SEGDUP_BEDPE"), sep="\t", header=None, usecols=range(1, 7),
                    names=["chr1", "start1", "end1", "chr2", "start2", "end2"])   # column 0 is a row index

# Members of the eight gene families: taken from merged_paralog_singlecopy_families.tsv (the
# same table behind Figure 5D), so C and D show the same genes. (The notebook matched gene
# names by substring, which also caught unrelated genes such as Casp10, Casp14, Ifngr1/2.)
fam_tbl = pd.read_csv(os.path.join(here, "merged_paralog_singlecopy_families.tsv"), sep="\t")
fam = (fam_tbl[fam_tbl["family"].isin(GENE_COLORS)]
       .rename(columns={"gene_name": "name"})[["chrom", "start", "end", "family", "name"]]
       .astype({"start": int, "end": int}))
assert set(fam["family"]) == set(GENE_COLORS), f"families missing: {set(GENE_COLORS) - set(fam['family'])}"
print(f"family genes: {len(fam)}")
print(fam.groupby("family", sort=False)["name"].apply(lambda x: ", ".join(x)).to_string())

# SegDup links overlapping any family gene; colour = family of the (first) gene they contain
links = []
for _, g in fam.iterrows():
    m = (((bedpe.chr1 == g.chrom) & (bedpe.start1 < g.end) & (bedpe.end1 > g.start)) |
         ((bedpe.chr2 == g.chrom) & (bedpe.start2 < g.end) & (bedpe.end2 > g.start)))
    for i in bedpe.index[m]:
        links.append((i, g.family))
link_family = {}
for i, f in links:                                   # keep first family per link, in GENE_COLORS order
    link_family.setdefault(i, f)
link_df = bedpe.loc[list(link_family)].assign(family=[link_family[i] for i in link_family])
link_df = link_df[link_df.chr1.isin(chr_sizes) & link_df.chr2.isin(chr_sizes)]
print(f"links drawn: {len(link_df)}; by family: {link_df.family.value_counts().to_dict()}")

# ---------------------------------------------------------------- plot
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "legend.fontsize": FONT_SIZE,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "savefig.bbox": None,      # pyCirclize sets "tight" on import, which resizes the panel
})
circos = Circos(sectors=chr_sizes, start=4, end=356, space=1.2, endspace=True)
for sector in circos.sectors:
    axis = sector.add_track((96, 99))
    axis.axis(fc="none", ec="black", lw=0.4)
    axis.xticks_by_interval(interval=50_000_000, tick_length=1.6, outer=True, show_label=False,
                            line_kws=dict(ec="black", lw=0.4))
    sector.text(sector.name.replace("chr", ""), r=106, size=FONT_SIZE, adjust_rotation=True,
                orientation="horizontal")

    dens = sector.add_track((84, 95))
    d = density[density["Chromosome"] == sector.name]
    if len(d):
        dens.bar(x=d["Start"].to_numpy(), height=d["percent_segdup"].to_numpy(), width=5_000_000,
                 align="edge", vmax=100.0, color=DENSITY_COLOR, lw=0)

    hg = sector.add_track((74, 83))
    rows = hot_genes[hot_genes["Chromosome"] == sector.name]
    hg.genomic_features([SeqFeature(FeatureLocation(int(r.Start), int(r.End))) for r in rows.itertuples()],
                        fc=HOTSPOT_GENE_COLOR, ec=HOTSPOT_GENE_COLOR, lw=0.3)

    gt = sector.add_track((64, 73))
    for g in fam[fam["chrom"] == sector.name].itertuples():
        # one radial tick per gene (0.8 pt): a gene is far narrower than a line at this scale
        mid = (g.start + g.end) / 2
        gt.line([mid, mid], [0, 1], vmin=0, vmax=1, color=GENE_COLORS[g.family], lw=0.8)

for r in link_df.itertuples():
    # One stroked Bezier curve between the two region midpoints (no fill). Filled pyCirclize
    # ribbons this thin are rendered by Illustrator as triangles.
    circos.link_line((r.chr1, (r.start1 + r.end1) / 2), (r.chr2, (r.start2 + r.end2) / 2), r1=63, r2=63,
                     color=GENE_COLORS[r.family], lw=0.6, alpha=0.9)

# Own figure + polar axes at an exact position, so the circle size is set in inches.
# pyCirclize draws r = 0-105 inside the axes; chromosome numbers (r = 106) sit just outside.
LABEL_PAD_IN = 0.16                                  # room for the chromosome numbers around the circle
KEY_IN = 0.36                                        # family key below
circle_in = PANEL_W - 2 * LABEL_PAD_IN
fig = plt.figure(figsize=(PANEL_W, PANEL_H))
ax = fig.add_axes((LABEL_PAD_IN / PANEL_W, (KEY_IN + LABEL_PAD_IN) / PANEL_H,
                   circle_in / PANEL_W, circle_in / PANEL_H), projection="polar")
circos.plotfig(ax=ax)

handles = [mpatches.Patch(color=c, label=g) for g, c in GENE_COLORS.items()]
leg = fig.legend(handles=handles, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.0), frameon=False,
                 handlelength=1.1, handleheight=0.8, columnspacing=1.1, handletextpad=0.4, labelspacing=0.3)
for t in leg.get_texts():
    t.set_fontstyle("italic")

fig.canvas.draw()
rend = fig.canvas.get_renderer()
circle_bottom = ax.get_window_extent(rend).y0
assert leg.get_window_extent(rend).y1 < circle_bottom + 0.06 * fig.dpi, "family key overlaps the circle"
for t in ax.texts:
    b = t.get_window_extent(rend)
    assert b.x0 >= 0 and b.x1 <= PANEL_W * fig.dpi and b.y1 <= PANEL_H * fig.dpi, "chromosome label off the panel"

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(f"{out_stem}.{{pdf,svg,png}}  ({opts.width_mm:.0f} x {opts.height_mm:.0f} mm)")
