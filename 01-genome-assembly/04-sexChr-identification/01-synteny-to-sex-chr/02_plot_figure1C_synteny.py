#!/usr/bin/env python
"""Figure 1C: mouse (GRCm39) vs degu macrosynteny, with mouse X highlighted.

Re-draws the JCVI karyotype plot at its printed size for Genome Research
(8-10 pt text, color key inside the panel, no titles). The synteny itself is not
recomputed: the script reads the seqids, BED and anchors files that
01_usefulMCscanScripts.sh produced, and uses JCVI's own ribbon drawing.
Chromosome/scaffold labels and the key are drawn here so their size is controlled.

Degu tracks are HapHiC groups; groupN is scaffold N of the assembly.

Run with the JCVI env:
    $ENV_JCVI/bin/python 02_plot_figure1C_synteny.py [--highlight-chr1]
"""

import argparse
import os
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import Bbox

import jcvi.graphics.karyotype as jk


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
LABEL_SIZE = 9          # track names ("Mouse", "Degu")
PANEL_WIDTHS = {"column": 3.31, "full": 6.89}   # in; GR one column 84 mm, full width 175 mm
PANEL_H = 1.97          # in; 50 mm
BAR_H = 0.03            # chromosome bar height, axes fraction (~1.5 mm)
RIBBON_LW = 0.25        # pt; edge of highlighted ribbons (JCVI default is 1 pt)

KEY_GAP_IN = 0.07       # in (~1.8 mm); scaffold labels to colour key
PAD_IN = 0.02           # in; margin kept above and below the content when trimming
LABEL_GAP_IN = 0.04     # in; min gap between neighbouring labels (~1 mm)
LABEL_OFFSET = 0.05    # label centre to bar centre, axes fraction
LABEL_OFFSET_LEADER = 0.09   # ... when the track has leader lines, to leave room for them
LEADER_MIN_IN = 0.02    # in; draw a leader line when a label moved more than this

COL_X = "#D55E00"       # Okabe-Ito vermillion: mouse X
COL_CHR1 = "#0072B2"    # Okabe-Ito blue: mouse Chr1 (only with --highlight-chr1)
COL_TRACK = "#BDBDBD"   # chromosome bars

# Prefixes written by highlight_chrom1.sh in the original run: r* = mouse ChrX, g* = mouse Chr1
SRC_PREFIX = {"r": "x", "g": "chr1"}

ap = argparse.ArgumentParser()
ap.add_argument("--highlight-chr1", action="store_true",
                help="also color the mouse Chr1 blocks (green in the submitted figure)")
ap.add_argument("--width", choices=PANEL_WIDTHS, default="column",
                help="column = 84 mm (labels in two rows), full = 175 mm (one row)")
opts = ap.parse_args()
PANEL_W = PANEL_WIDTHS[opts.width]
X0 = 0.46 / PANEL_W         # tracks start 0.46 in from the left, leaving room for "Mouse"/"Degu"
X1 = 1 - 0.08 / PANEL_W     # and end 0.08 in from the right
LABEL_GAP = LABEL_GAP_IN / PANEL_W

in_dir = _cfg("JCVI_KARYOTYPE_DIR")
out_dir = _cfg("FIG1C_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
tag = "chrX_chr1" if opts.highlight_chr1 else "chrX"
out_stem = os.path.join(out_dir, f"figure1C_mouse_degu_synteny_{tag}_{opts.width}")

# --- inputs for JCVI, written to out_dir so the original run stays untouched ---
colors = {"x": COL_X, "chr1": COL_CHR1 if opts.highlight_chr1 else None}
simple = os.path.join(out_dir, f"mouse.degus.anchors.{tag}.simple")
with open(os.path.join(in_dir, "mouse.degus.anchors.chr1.chrX.simple1")) as fin, open(simple, "w") as fout:
    for row in fin:
        if "*" in row:
            prefix, row = row.split("*", 1)
            col = colors[SRC_PREFIX[prefix]]
            if col:
                row = f"{col}*{row}"
        fout.write(row)

layout = os.path.join(out_dir, f"layout.{tag}")
with open(layout, "w") as f:
    f.write("# y, xstart, xend, rotation, color, label, va, bed\n")
    f.write(f" 0.80, {X0:.4f}, {X1:.4f}, 0, {COL_TRACK}, Mouse, top, {in_dir}/mouse.bed\n")
    f.write(f" 0.30, {X0:.4f}, {X1:.4f}, 0, {COL_TRACK}, Degu, bottom, {in_dir}/degus.bed\n")
    f.write("# edges\n")
    f.write(f"e, 0, 1, {simple}\n")

# --- plot ---
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})



class CleanBar:
    """Stand-in for JCVI's HorizontalChromosome: one filled, rounded bar, no outline."""

    def __init__(self, ax, x1, x2, y, height=BAR_H, fc=None, **_):
        x1, x2 = sorted((x1, x2))
        aspect = PANEL_W / PANEL_H             # makes the rounded ends circular on paper
        r = min(height / 2 / aspect, (x2 - x1) / 2)
        self.patch = FancyBboxPatch((x1, y - height / 2), x2 - x1, height,
                                    boxstyle=f"round,pad=0,rounding_size={r}",
                                    mutation_aspect=aspect, fc=fc or COL_TRACK, lw=0, zorder=3)
        ax.add_patch(self.patch)

    def set_transform(self, tr):
        self.patch.set_transform(tr)


def thin_shade(*args, **kwargs):
    if kwargs.get("lw"):
        kwargs["lw"] = RIBBON_LW
    return _Shade(*args, **kwargs)


_Shade = jk.Shade
jk.HorizontalChromosome = CleanBar
jk.Shade = thin_shade

fig = plt.figure(figsize=(PANEL_W, PANEL_H))
root = fig.add_axes((0, 0, 1, 1))
k = jk.Karyotype(root, os.path.join(in_dir, "seqids"), layout,
              height=BAR_H, plot_label=False, plot_circles=False)


def short(sid):
    # "Chr7" -> "7", "group12" -> "12"
    return sid.replace("Chr", "").replace("group", "")


# Labels, horizontal (rotating to 45 deg saves no width: parallel slanted labels need
# (text height + gap) / sin 45 of horizontal spacing, about the width of a 2-digit label).
# Labels stay in one row; any label pushed sideways to make room gets a leader line.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()


def spread(xs, w, right=0.998):
    """Push labels apart so neighbours keep LABEL_GAP, then pull the run back from the right edge."""
    xs = list(xs)
    for i in range(1, len(xs)):
        xs[i] = max(xs[i], xs[i - 1] + (w[i - 1] + w[i]) / 2 + LABEL_GAP)
    for i in range(len(xs) - 1, -1, -1):
        limit = right - w[i] / 2 if i == len(xs) - 1 else xs[i + 1] - (w[i] + w[i + 1]) / 2 - LABEL_GAP
        xs[i] = min(xs[i], limit)
    return xs


# Only chromosomes/scaffolds that carry a highlighted (colored) synteny block are labeled.
involved = [set(), set()]
for i, j, blocks, _ in k.layout.edges:
    for a, b, c, d, _, _, hl in blocks:
        if hl:
            involved[i] |= {k.tracks[i].order_in_chr[g][0] for g in (a, b) if g in k.tracks[i].order_in_chr}
            involved[j] |= {k.tracks[j].order_in_chr[g][0] for g in (c, d) if g in k.tracks[j].order_in_chr}

for tr, sign, keep in zip(k.tracks, (+1, -1), involved):
    texts, x0s, ws = [], [], []
    for sid in (sid for sid in tr.seqids if sid in keep):
        x = tr.offsets[sid] + tr.ratio * tr.sizes[sid] / 2
        t = root.text(x, tr.y, short(sid), ha="center", va="center",
                      fontsize=FONT_SIZE, transform=root.transAxes)
        texts.append(t)
        x0s.append(x)
        ws.append(t.get_window_extent(renderer).transformed(root.transAxes.inverted()).width)
    # Always one row: labels pushed sideways get a leader line to their chromosome.
    xs = spread(x0s, ws)
    moved = [abs(x - x0) * PANEL_W > LEADER_MIN_IN for x, x0 in zip(xs, x0s)]
    label_y = tr.y + sign * (LABEL_OFFSET_LEADER if any(moved) else LABEL_OFFSET)
    text_half_h = FONT_SIZE * 0.36 / 72 / PANEL_H        # half a digit's height, axes fraction
    for t, x, x0, m in zip(texts, xs, x0s, moved):
        t.set_position((x, label_y))
        if m:
            root.plot([x0, x], [tr.y + sign * (BAR_H / 2 + 0.006), label_y - sign * (text_half_h + 0.008)],
                      color="#4D4D4D", lw=0.5, solid_capstyle="butt", transform=root.transAxes)

    root.text(0.005, tr.y, tr.label, ha="left", va="center",
              fontsize=LABEL_SIZE, transform=root.transAxes)

# Color key, inside the panel, KEY_GAP_IN below the lowest scaffold label or leader line
fig.canvas.draw()
to_ax = root.transAxes.inverted()
lowest = min(to_ax.transform(a.get_window_extent(renderer))[:, 1].min()
             for a in root.texts + root.lines if a.get_visible())
handles = [Line2D([], [], color=COL_X, lw=3, label="Mouse X")]
if opts.highlight_chr1:
    handles.append(Line2D([], [], color=COL_CHR1, lw=3, label="Mouse 1"))
handles.append(Line2D([], [], color="gainsboro", lw=3, label="Other"))
root.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.55, lowest - KEY_GAP_IN / PANEL_H),
            ncol=len(handles), frameon=False, fontsize=FONT_SIZE, handlelength=1.4,
            columnspacing=1.2, borderaxespad=0, borderpad=0)

root.set_xlim(0, 1)
root.set_ylim(0, 1)
root.set_axis_off()

# Trim empty space above and below; the width stays exactly PANEL_W.
fig.canvas.draw()
drawn = Bbox.union([a.get_window_extent(renderer) for a in root.texts + root.lines + root.patches
                    + [root.get_legend()]]).transformed(fig.dpi_scale_trans.inverted())   # inches
crop = Bbox([[0, max(drawn.y0 - PAD_IN, 0)], [PANEL_W, min(drawn.y1 + PAD_IN, PANEL_H)]])
print(f"panel: {PANEL_W * 25.4:.0f} x {crop.height * 25.4:.1f} mm ({PANEL_W} x {crop.height:.2f} in)")

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600, bbox_inches=crop)
print(out_stem + ".{pdf,svg,png}")
