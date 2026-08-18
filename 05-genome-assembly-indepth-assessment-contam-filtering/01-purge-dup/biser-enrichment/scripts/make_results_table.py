#!/usr/bin/env python3
"""Generate publication-ready enrichment results table matching the style of
category_exploration_v4_busco.ipynb Panel J."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties
import json, re
from pathlib import Path

OUTPUT_ROOT = Path(__file__).resolve().parent.parent
TABLES_DIR = OUTPUT_ROOT / "tables"
FIGURES_DIR = OUTPUT_ROOT / "figures"

# ── Colors (from purge_dup_distribution.ipynb) ─────────────────────────
CAT_COLORS = {
    "REPEAT":   "#FDAE61",
    "HAPLOTIG": "#D73027",
    "JUNK":     "#4D4D4D",
    "OVLP":     "#1B9E77",
    "HIGHCOV":  "#762A83",
}
CAT_ORDER = ["HAPLOTIG", "HIGHCOV", "JUNK", "OVLP", "REPEAT"]

# ── Load data ──────────────────────────────────────────────────────────
def load_table(path):
    with open(path) as f:
        headers = f.readline().strip().split("\t")
        return [dict(zip(headers, l.strip().split("\t"))) for l in f]

obs = {r["category"]: r for r in load_table(TABLES_DIR / "purge_dups_biser_observed.tsv")}
perm = {r["category"]: r for r in load_table(TABLES_DIR / "purge_dups_biser_permutation_enrichment.tsv")}
n_perm = perm[CAT_ORDER[0]]["n_permutations"]

# ── Build table rows ───────────────────────────────────────────────────
rows = []
for cat in CAT_ORDER:
    o, p = obs[cat], perm[cat]
    fold = float(p["fold_enrichment"])
    pval = float(p["empirical_p_enrichment"])
    fdr = float(p["fdr_bh"])
    sig = "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*" if fdr < 0.05 else "ns"
    rows.append({
        "cat": cat,
        "color": CAT_COLORS[cat],
        "cells": [
            cat,
            o["n_intervals"],
            f'{float(o["total_mb"]):.1f}',
            f'{float(o["pct_bases_overlapping"]):.1f}%',
            f'{float(p["observed_overlap_mb"]):.2f}',
            f'{float(p["null_mean_overlap_bp"])/1e6:.2f}',
            f'{fold:.1f}×',
            f'{pval:.4f}',
            f'{fdr:.4f}',
            sig,
        ],
    })

col_labels = [
    "Category",   "N intervals", "Total\n(Mb)", "Raw BISER\n(%)",
    "Observed\n(Mb)", "Null mean\n(Mb)", "Fold\nenrichment", "P (enr)",
    "FDR", "Sig.",
]
# Narrow count columns, wider label columns
col_widths = [0.14, 0.06, 0.06, 0.09, 0.09, 0.10, 0.10, 0.10, 0.10, 0.06,
              0.10]  # 11 values = sum to 1.0
# Actually let me simplify: 10 columns
col_labels10 = [
    "Category", "N", "Total (Mb)", "Raw BISER (%)",
    "Observed (Mb)", "Null mean (Mb)", "Fold enrichment", "P (enr)",
    "FDR", "Sig.",
]
col_widths10 = [0.145, 0.055, 0.075, 0.095, 0.10, 0.105, 0.095, 0.105, 0.105, 0.12]
assert abs(sum(col_widths10) - 1.0) < 1e-6

# ── Typography ─────────────────────────────────────────────────────────
FIG_W = 14.0
TITLE_FS, HEADER_FS, CAT_FS, BODY_FS, COUNT_FS = 18, 13, 11.5, 10.5, 11
HEADER_H_PT, ROW_PAD_PT = 40, 5
TOP_MARGIN, TITLE_H, TITLE_GAP, BOT_MARGIN = 6, 24, 6, 5

body_font   = FontProperties(family="sans-serif", size=BODY_FS)
cat_font    = FontProperties(family="sans-serif", size=CAT_FS, weight="bold")
header_font = FontProperties(family="sans-serif", size=HEADER_FS, weight="bold")
title_font  = FontProperties(family="sans-serif", size=TITLE_FS, weight="bold")
count_font  = FontProperties(family="sans-serif", size=COUNT_FS, weight="bold")

def norm(s):
    return re.sub(r"\s+", " ", str(s)).strip()

def text_width(renderer, text, fp):
    w, _, _ = renderer.get_text_width_height_descent(text, fp, ismath=False)
    return w

def wrap(text, max_px, renderer, fp):
    text = norm(text)
    if not text: return ""
    words = text.split()
    lines, cur = [], words[0]
    for w in words[1:]:
        cand = f"{cur} {w}"
        if text_width(renderer, cand, fp) <= max_px:
            cur = cand
        else:
            lines.append(cur); cur = w
    lines.append(cur)
    return "\n".join(lines)

# ── Measure & wrap ─────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIG_W, 6), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
fig.canvas.draw()
renderer = fig.canvas.get_renderer()

fig_w_px = FIG_W * fig.dpi
table_w_px = fig_w_px * 0.96
h_pad_px = 6.5 / 72 * fig.dpi

wrapped_rows = []
row_lines = []
for row in rows:
    wrapped = []
    for ci, val in enumerate(row["cells"]):
        cw = table_w_px * col_widths10[ci]
        avail = cw - 2 * h_pad_px - 5
        fp = cat_font if ci == 0 else count_font if ci == 1 else body_font
        wrapped.append(wrap(val, avail, renderer, fp))
    wrapped_rows.append({"text": wrapped, "color": row["color"]})
    row_lines.append(max(t.count("\n") + 1 for t in wrapped))

line_h_pt = BODY_FS * 1.15
row_heights_pt = [l * line_h_pt + 2 * ROW_PAD_PT for l in row_lines]
table_h_pt = HEADER_H_PT + sum(row_heights_pt)
fig_h_pt = TOP_MARGIN + TITLE_H + TITLE_GAP + table_h_pt + BOT_MARGIN
fig_h_in = fig_h_pt / 72
fig.set_size_inches(FIG_W, fig_h_in, forward=True)
fig.canvas.draw()

def pt2y(pt):
    return pt / fig_h_pt

title_y = 1.0 - pt2y(TOP_MARGIN)
table_top = 1.0 - pt2y(TOP_MARGIN + TITLE_H + TITLE_GAP)
header_h = pt2y(HEADER_H_PT)
row_hs = [pt2y(h) for h in row_heights_pt]
h_pad_f = (6.5 / 72) / FIG_W

# Column edges
col_lefts = [0.02]
for i in range(len(col_widths10) - 1):
    col_lefts.append(col_lefts[-1] + 0.96 * col_widths10[i])
col_disp_w = [0.96 * w for w in col_widths10]

# ── Title ──────────────────────────────────────────────────────────────
ax.text(0.5, title_y,
    f"purge_dups category enrichment for BISER segmental duplications\n"
    f"({int(n_perm):,} chromosome- and length-matched permutations, BH-FDR corrected)",
    fontproperties=title_font, color="black", ha="center", va="top")

# ── Header ─────────────────────────────────────────────────────────────
hb = table_top - header_h
for ci, lbl in enumerate(col_labels10):
    r = Rectangle((col_lefts[ci], hb), col_disp_w[ci], header_h,
                   facecolor="#333333", edgecolor="#666666", linewidth=0.7)
    ax.add_patch(r)
    ax.text(col_lefts[ci] + col_disp_w[ci] / 2, hb + header_h / 2, lbl,
            fontproperties=header_font, color="white", ha="center", va="center")

# ── Body rows ──────────────────────────────────────────────────────────
cur_top = hb
for ri, row in enumerate(wrapped_rows):
    rb = cur_top - row_hs[ri]
    for ci, txt in enumerate(row["text"]):
        cl, cw = col_lefts[ci], col_disp_w[ci]
        if ci == 0:
            fc, tc, fp, ha, tx = row["color"], "white", cat_font, "left", cl + h_pad_f
        elif ci == 1:
            fc, tc, fp, ha, tx = "#FAFAFA", "#222222", count_font, "center", cl + cw / 2
        elif ci == 9:
            fc, tc, fp, ha, tx = "#FAFAFA", "#222222", count_font, "center", cl + cw / 2
        else:
            fc, tc, fp, ha, tx = "#FAFAFA", "#222222", body_font, "left", cl + h_pad_f

        r = Rectangle((cl, rb), cw, row_hs[ri], facecolor=fc,
                       edgecolor="#666666", linewidth=0.7)
        ax.add_patch(r)
        ax.text(tx, rb + row_hs[ri] / 2, txt, fontproperties=fp,
                color=tc, ha=ha, va="center", linespacing=1.05, clip_on=True)
    cur_top = rb

out = FIGURES_DIR / "enrichment_results_table.png"
fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.03)
plt.close(fig)
print(f"Saved: {out}")
