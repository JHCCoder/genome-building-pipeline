#!/usr/bin/env python3
"""
Replace the OctDeg2.0 rows of the Figure 3D tables with the values from the matched-annotation
Cell Ranger ARC reruns (07_cellranger_arc_count_octDeg2.sh, Liftoff-only reference), for the
samples in data/fig3d_multiome_libraries.tsv. OctDeg1.0 rows are left unchanged (they come from
the original OctDeg1 runs listed in that table's octdeg1_run column).

Each table is backed up to <name>.pre-<YYYYMMDD> before it is rewritten.

Inputs: $FIG3D_CR_DIR/count/<sample>_OctDeg2_liftoffOnly/outs/summary.csv
Tables (data/): feature_linkage.txt, gene_per_cell.txt, fragment_per_cell.txt

Run with the python-visualizations env:
    python 08_update_fig3d_octdeg2_metrics.py
"""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd


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


DATA = Path(__file__).resolve().parent / "data"
COUNT_DIR = Path(_cfg("FIG3D_CR_DIR")) / "count"
SAMPLES = pd.read_csv(DATA / "fig3d_multiome_libraries.tsv", sep="\t")["sample"].tolist()

# table -> (value column in the table, column in Cell Ranger ARC summary.csv)
TABLES = {
    "feature_linkage.txt": ("Feature linkages detected", "Feature linkages detected"),
    "gene_per_cell.txt": ("Median genes per cell", "GEX Median genes per cell"),
    "fragment_per_cell.txt": ("Median high quality fragment per cell",
                              "ATAC Median high-quality fragments per cell"),
}

summaries = {}
for s in SAMPLES:
    f = COUNT_DIR / f"{s}_OctDeg2_liftoffOnly" / "outs" / "summary.csv"
    if not f.is_file():
        sys.exit(f"ERROR: missing {f}")
    row = pd.read_csv(f).iloc[0]
    if row["Genome"] != "OctDeg2_liftoffOnly":
        sys.exit(f"ERROR: {f} was run on genome {row['Genome']}, expected OctDeg2_liftoffOnly")
    summaries[s] = row

stamp = time.strftime("%Y%m%d")
for fname, (col, summary_col) in TABLES.items():
    path = DATA / fname
    tbl = pd.read_csv(path, sep="\t")
    backup = path.with_name(f"{path.name}.pre-{stamp}")
    if not backup.exists():
        shutil.copy2(path, backup)
    print(f"== {fname}  (backup: {backup.name})")
    for s in SAMPLES:
        new = int(round(float(summaries[s][summary_col])))
        mask = (tbl["samples"] == s) & (tbl["genome"] == "octDeg2")
        if not mask.any():
            sys.exit(f"ERROR: no octDeg2 row for {s} in {fname}")
        old = tbl.loc[mask, col].unique().tolist()
        tbl.loc[mask, col] = new
        print(f"   {s:10s} octDeg2  {old} -> {new}")
    tbl.to_csv(path, sep="\t", index=False)
