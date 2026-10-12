#!/usr/bin/env python
"""Figure 2C: five-category breakdown of the protein-coding genes in the final annotation.

Category logic is copied from 01_annotation_pie_chart.ipynb (cells 4-7):
  LOC genes                          gene ID matches LOC + 9 digits
  paralog copy                       gene ID ends in -L1 / -DL1 / -RL1 ...
  novel                              suffix-stripped name is absent from the
                                     protein-coding genes transferred by Liftoff
Drawn at printed size for Genome Research (8 pt text, key inside the panel, no title).
The counts are also written to figure2C_counts.tsv.

Run with the python-visualizations env:
    python 02_plot_figure2C_annotation_pie.py [--width-mm 55]
"""

import argparse
import os
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox


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


FONT_SIZE = 8           # pt at final print size; Genome Research: 8-10 pt
PANEL_H = 3.0           # in; working canvas, trimmed to the content when saved
PIE_D_IN = 1.30         # pie diameter, in (~33 mm)
PARALOG = r"-(?:D?R?L\d+)$"

# (category, key text, face color, hatch). Colours from the original figure
# (01_annotation_pie_chart.ipynb): light blue = characterized, light red = novel,
# light green = LOC; paralog copies share their parent's colour with a dot pattern.
CATEGORIES = [
    ("Characterized genes", "Characterized genes", "#aec7e8", None),
    ("Characterized-gene paralog copies", "Characterized-gene\nparalog copies", "#aec7e8", "...."),
    ("Novel genes", "Novel genes", "#ff9896", None),
    ("Novel-gene paralog copies", "Novel-gene\nparalog copies", "#ff9896", "...."),
    ("LOC genes", "LOC genes", "#98df8a", None),
]
EDGE = "black"

ap = argparse.ArgumentParser()
ap.add_argument("--width-mm", type=float, default=55)
opts = ap.parse_args()
PANEL_W = opts.width_mm / 25.4

out_dir = _cfg("FIG2C_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, f"figure2C_annotation_pie_{opts.width_mm:.0f}mm")

# --- counts (notebook logic) ---
cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
lo = pd.read_csv(_cfg("FIG2C_LIFTOFF_GFF"), sep="\t", comment="#", header=None, names=cols,
                 usecols=["type", "attributes"])
lo_biotype = lo["attributes"].str.extract(r"gene_biotype=([^;]+)", expand=False)
lo_name = lo["attributes"].str.extract(r"ID=gene-([^;]+)", expand=False)
reference = set(lo_name[lo_biotype == "protein_coding"].dropna().str.upper())

final = pd.read_csv(_cfg("FINAL_ANNOTATION_GFF"), sep="\t", comment="#", header=None, names=cols,
                    usecols=["type", "attributes"])
genes = final[final["type"] == "gene"].copy()
genes["biotype"] = genes["attributes"].str.extract(r"gene_biotype=([^;]+)", expand=False).fillna("protein_coding")
genes["gene_name"] = genes["attributes"].str.extract(r"ID=gene-([^;]+)", expand=False).str.upper()
genes["is_loc"] = genes["gene_name"].str.match(r"^LOC\d{9}$", na=False)
genes["novel"] = ~genes["gene_name"].str.replace(PARALOG, "", regex=True).isin(reference)
genes = genes[["biotype", "gene_name", "is_loc", "novel"]].dropna().drop_duplicates()
genes = genes[genes["biotype"] == "protein_coding"]
genes["paralog"] = genes["gene_name"].str.contains(PARALOG, na=False, regex=True)
genes["category"] = np.select(
    [genes["is_loc"], genes["paralog"] & genes["novel"], genes["paralog"], genes["novel"]],
    ["LOC genes", "Novel-gene paralog copies", "Characterized-gene paralog copies", "Novel genes"],
    default="Characterized genes")
counts = genes["category"].value_counts().reindex([c[0] for c in CATEGORIES]).fillna(0).astype(int)
total = counts.sum()
pd.DataFrame({"category": counts.index, "count": counts.values,
              "percent": (counts.values / total * 100).round(1)}).to_csv(
    os.path.join(out_dir, "figure2C_counts.tsv"), sep="\t", index=False)
print(counts.to_string(), f"\ntotal {total:,}")

# --- plot ---
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE, "hatch.linewidth": 0.5, "hatch.color": EDGE,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
fig = plt.figure(figsize=(PANEL_W, PANEL_H))

# Pie at the top, centred
pie_ax = fig.add_axes(((PANEL_W - PIE_D_IN) / 2 / PANEL_W, 1 - (PIE_D_IN + 0.05) / PANEL_H,
                       PIE_D_IN / PANEL_W, PIE_D_IN / PANEL_H))
wedges, _ = pie_ax.pie(counts.values, startangle=90, counterclock=False,
                       colors=[c[2] for c in CATEGORIES],
                       explode=[0.06 if "Novel" in c[0] else 0 for c in CATEGORIES],
                       wedgeprops=dict(lw=0.5, edgecolor=EDGE))
for w, (_, _, _, hatch) in zip(wedges, CATEGORIES):
    if hatch:
        w.set_hatch(hatch)
pie_ax.set_aspect("equal")

# Key under the pie: swatch | category | count | percent
key_ax = fig.add_axes((0, 0, 1, 1 - (PIE_D_IN + 0.0) / PANEL_H))
key_ax.set_xlim(0, PANEL_W)
key_ax.set_ylim(0, 1)
key_ax.set_axis_off()
line_in = FONT_SIZE * 1.2 / 72
key_h_in = (1 - (PIE_D_IN + 0.0) / PANEL_H) * PANEL_H
y_in = key_h_in - 0.02
SW = 0.09                                   # swatch size, in
x_sw, x_name, x_n, x_pct = 0.03, 0.17, PANEL_W - 0.42, PANEL_W - 0.03
for (cat, text, color, hatch), n in zip(CATEGORIES, counts.values):
    nlines = text.count("\n") + 1
    yc = y_in - line_in * nlines / 2
    key_ax.add_patch(Rectangle((x_sw, (yc - SW / 2) / key_h_in), SW, SW / key_h_in,
                               facecolor=color, hatch=hatch, edgecolor=EDGE, lw=0.5))
    key_ax.text(x_name, yc / key_h_in, text, ha="left", va="center", linespacing=1.0)
    key_ax.text(x_n, yc / key_h_in, f"{n:,}", ha="right", va="center")
    key_ax.text(x_pct, yc / key_h_in, f"{n / total * 100:.1f}%", ha="right", va="center")
    y_in -= line_in * nlines + 0.035

# Trim empty space above and below; the width stays exactly PANEL_W.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
drawn = Bbox.union([a.get_window_extent(renderer) for a in wedges + key_ax.texts + key_ax.patches]
                   ).transformed(fig.dpi_scale_trans.inverted())
crop = Bbox([[0, drawn.y0 - 0.02], [PANEL_W, drawn.y1 + 0.02]])
print(f"panel: {opts.width_mm:.0f} x {crop.height * 25.4:.1f} mm ({PANEL_W:.2f} x {crop.height:.2f} in)")

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600, bbox_inches=crop)
print(out_stem + ".{pdf,svg,png}")
