#!/usr/bin/env python
"""Figure 1F: BLAST hit of grey squirrel Sry (XM_047537220.1) on degu scaffold 23.

Reads the BLAST XML written by 04_blastSryGenes.sh and draws the query (top) and
scaffold 23 (bottom, reversed because the hit is on the minus strand) with the
aligned region joined, plus the hit statistics as live text. Replaces the earlier
Kablammo-style panel, whose statistics were a ~200 dpi screenshot.

Run with the python-visualizations env:
    python 04_plot_figure1F_sry_blast.py [--width column|full]
"""

import argparse
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


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
PANEL_WIDTHS = {"column": 3.31, "full": 6.89}   # in; GR one column 84 mm, full width 175 mm
PANEL_H = 1.45          # in (~37 mm)
TRACKS_FRAC = {"column": 0.44, "full": 0.70}    # share of the width used by the two tracks
QUERY_ID = "XM_047537220.1"                     # grey squirrel (Sciurus carolinensis) Sry
SCAFFOLD_LABEL = "23"
HIT_COLOR = "#0072B2"   # Okabe-Ito blue

ap = argparse.ArgumentParser()
ap.add_argument("--width", choices=PANEL_WIDTHS, default="column")
opts = ap.parse_args()
PANEL_W = PANEL_WIDTHS[opts.width]

out_dir = _cfg("FIG1F_OUT_DIR")
os.makedirs(out_dir, exist_ok=True)
out_stem = os.path.join(out_dir, f"figure1F_sry_blast_{opts.width}")

# --- the best HSP for the query ---
root = ET.parse(_cfg("SRY_BLAST_XML")).getroot()
hsp = None
for it in root.iter("Iteration"):
    if it.findtext("Iteration_query-def").startswith(QUERY_ID):
        qlen = int(it.findtext("Iteration_query-len"))
        hit = it.find("Iteration_hits/Hit")
        slen = int(hit.findtext("Hit_len"))
        hsp = max(hit.iter("Hsp"), key=lambda h: float(h.findtext("Hsp_bit-score")))
assert hsp is not None, f"{QUERY_ID} has no hit in the BLAST XML"
v = {k: hsp.findtext(f"Hsp_{k}") for k in
     ("bit-score", "evalue", "query-from", "query-to", "hit-from", "hit-to", "hit-frame",
      "align-len", "identity")}
q0, q1 = int(v["query-from"]), int(v["query-to"])
s0, s1 = sorted((int(v["hit-from"]), int(v["hit-to"])))
minus = int(v["hit-frame"]) < 0

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans"],
    "font.size": FONT_SIZE,
    "mathtext.fontset": "custom", "mathtext.rm": "Nimbus Sans", "mathtext.it": "Nimbus Sans:italic", "mathtext.bf": "Nimbus Sans:bold", "mathtext.cal": "Nimbus Sans", "mathtext.sf": "Nimbus Sans",
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

fig = plt.figure(figsize=(PANEL_W, PANEL_H))
ax = fig.add_axes((0, 0, 1, 1))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_axis_off()

X0 = 0.03                                 # tracks start
X1 = X0 + TRACKS_FRAC[opts.width]         # tracks end
YQ, YS = 0.68, 0.30                       # query / scaffold track heights (axes fraction)
TICK = 0.035                              # tick length, axes fraction
LW = 0.75


def xq(bp):                               # query position -> x
    return X0 + (X1 - X0) * bp / qlen


def xs(bp):                               # scaffold position -> x (reversed for a minus-strand hit)
    f = bp / slen
    return X0 + (X1 - X0) * ((1 - f) if minus else f)


def track(y, ticks, to_x, labels, up):
    ax.plot([X0, X1], [y, y], color="black", lw=LW, solid_capstyle="butt")
    d = TICK if up else -TICK
    for t, lab in zip(ticks, labels):
        x = to_x(t)
        ax.plot([x, x], [y, y + d], color="black", lw=0.5)
        ax.text(x, y + d * 1.6, lab, ha="center", va="bottom" if up else "top", fontsize=FONT_SIZE)


track(YQ, range(0, qlen + 1, 200), xq, [f"{t}" for t in range(0, qlen + 1, 200)], up=True)
mb = range(0, int(slen / 1e6) + 1, 20)
track(YS, [m * 1e6 for m in mb], xs, [f"{m}" for m in mb], up=False)

# Units at the right end of each track
units = [ax.text(X1 + 0.012, YQ, "bp", ha="left", va="center", fontsize=FONT_SIZE),
         ax.text(X1 + 0.012, YS, "Mb", ha="left", va="center", fontsize=FONT_SIZE)]

# Track names
ax.text(X0, 0.98, r"Grey squirrel $\it{Sry}$ mRNA", ha="left", va="top", fontsize=FONT_SIZE)
ax.text(X0, 0.02, f"Degu scaffold {SCAFFOLD_LABEL}: {s0:,}\u2013{s1:,} ({'\u2212' if minus else '+'})",
        ha="left", va="bottom", fontsize=FONT_SIZE)

# Aligned region
ax.add_patch(Polygon([(xq(q0), YQ), (xq(q1), YQ), (xs(s1), YS), (xs(s0), YS)],
                     closed=True, fc=HIT_COLOR, ec=HIT_COLOR, lw=0.5, zorder=0))

# Hit statistics, live text
mant, exp = f"{float(v['evalue']):.1e}".split("e")
ident = int(v["identity"]) / int(v["align-len"]) * 100
stats = [
    ("Bit score", f"{float(v['bit-score']):.1f}"),
    ("E-value", f"{mant}e\u2212{abs(int(exp))}"),   # no superscript: it would print below 8 pt
    ("Identity", f"{v['identity']}/{v['align-len']} ({ident:.0f}%)"),
    ("Alignment", f"{v['align-len']} bp"),
    ("Query", f"{q0}\u2013{q1}"),
]
# Keys in bold, values aligned 1.5 mm right of the widest key; block centred between the tracks.
renderer = fig.canvas.get_renderer()
units_right = max(t.get_window_extent(renderer).x1 for t in units) / fig.dpi / PANEL_W
sx = units_right + 0.16 / PANEL_W         # 4 mm right of the "bp"/"Mb" units
line = FONT_SIZE * 1.3 / 72 / PANEL_H
y0 = (YQ + YS) / 2 + line * (len(stats) - 1) / 2
keys = [ax.text(sx, y0 - i * line, k, ha="left", va="center", fontsize=FONT_SIZE, fontweight="bold")
        for i, (k, _) in enumerate(stats)]
key_w = max(t.get_window_extent(renderer).width for t in keys) / fig.dpi / PANEL_W   # axes fraction
vx = sx + key_w + 0.06 / PANEL_W                                                        # + 1.5 mm
vals = [ax.text(vx, y0 - i * line, val, ha="left", va="center", fontsize=FONT_SIZE)
        for i, (_, val) in enumerate(stats)]
right = max(t.get_window_extent(renderer).x1 for t in vals) / fig.dpi
assert right <= PANEL_W, f"statistics run {right - PANEL_W:.2f} in past the panel edge"

for ext in ("pdf", "svg", "png"):
    path = f"{out_stem}.{ext}"
    if os.path.exists(path):
        os.replace(path, f"{path}.pre-{time.strftime('%Y%m%d-%H%M%S')}")
    fig.savefig(path, dpi=600)
print(out_stem + ".{pdf,svg,png}")
