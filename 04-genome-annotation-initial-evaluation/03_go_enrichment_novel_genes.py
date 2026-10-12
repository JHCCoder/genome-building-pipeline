#!/usr/bin/env python
"""GO enrichment of the novel genes against a degu background (input to Figure 2D/E).

Script version of 03_go_enrichment_novel_genes.ipynb, rebuilt on the current final
annotation ($FINAL_ANNOTATION_GFF) so the gene sets match Figure 2C.

Novel genes (same rule as Figure 2C / 02_plot_figure2C_annotation_pie.py):
  protein-coding genes of the final annotation whose suffix-stripped name
  (-L1 / -DL1 / -RL1 removed) is absent from the protein-coding genes transferred by
  Liftoff, excluding LOC genes.
Foreground: the unique gene symbols of the novel genes ("Novel genes" plus the parent
  symbols of "Novel-gene paralog copies"; a paralog copy has the same symbol as its parent).
Background: every gene symbol (suffix-stripped) in the final annotation.
Libraries: local Enrichr 2025 GO BP/MF/CC GMT files, so gseapy runs offline.

Writes $FIG2_NOVEL_GENES_TSV (one row per novel gene) and $FIG2_GO_TSV (all terms).

Run with the go-analysis env:
    python 03_go_enrichment_novel_genes.py
"""

import os
import re
import subprocess
import sys

import gseapy as gp
import numpy as np
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


PARALOG = r"-(?:D?R?L\d+)$"
LIBRARIES = ["human_GO_bp_2025.gmt", "human_GO_mf_2025.gmt", "human_GO_cc_2025.gmt"]

for out_var in ("FIG2_NOVEL_GENES_TSV", "FIG2_GO_TSV"):
    if os.path.exists(_cfg(out_var)):
        sys.exit(f"ERROR: {_cfg(out_var)} exists; move it aside (e.g. *.pre-<date>) before regenerating")

cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
lo = pd.read_csv(_cfg("FIG2C_LIFTOFF_GFF"), sep="\t", comment="#", header=None, names=cols,
                 usecols=["type", "attributes"])
reference = set(lo["attributes"].str.extract(r"ID=gene-([^;]+)", expand=False)
                [lo["attributes"].str.contains("gene_biotype=protein_coding", regex=False)].dropna().str.upper())

final = pd.read_csv(_cfg("FINAL_ANNOTATION_GFF"), sep="\t", comment="#", header=None, names=cols,
                    usecols=["type", "attributes"])
genes = final[final["type"] == "gene"].copy()
genes["biotype"] = genes["attributes"].str.extract(r"gene_biotype=([^;]+)", expand=False).fillna("protein_coding")
genes["gene_name"] = genes["attributes"].str.extract(r"ID=gene-([^;]+)", expand=False).str.upper()
genes = genes.dropna(subset=["gene_name"])
genes["symbol"] = genes["gene_name"].str.replace(PARALOG, "", regex=True)
genes["is_loc"] = genes["gene_name"].str.match(r"^LOC\d{9}$")
genes["paralog_copy"] = genes["gene_name"].str.contains(PARALOG, regex=True)
genes["novel"] = ~genes["symbol"].isin(reference)

background = sorted(set(genes["symbol"]))
pc = genes[genes["biotype"] == "protein_coding"][["gene_name", "symbol", "is_loc", "paralog_copy", "novel"]]
pc = pc.drop_duplicates()
novel = pc[pc["novel"] & ~pc["is_loc"]].copy()
novel["category"] = np.where(novel["paralog_copy"], "Novel-gene paralog copies", "Novel genes")
foreground = sorted(set(novel["symbol"]))

print(f"background symbols: {len(background):,}")
print(novel["category"].value_counts().to_string())
print(f"foreground (unique novel symbols): {len(foreground):,}")

os.makedirs(os.path.dirname(_cfg("FIG2_NOVEL_GENES_TSV")), exist_ok=True)
novel[["gene_name", "symbol", "category"]].sort_values("gene_name").to_csv(
    _cfg("FIG2_NOVEL_GENES_TSV"), sep="\t", index=False)

gmt_dir = _cfg("FIG2_GMT_DIR")
res = gp.enrichr(gene_list=foreground, gene_sets=[os.path.join(gmt_dir, g) for g in LIBRARIES],
                 organism="human", background=background, outdir=None).results
res["Gene_set"] = res["Gene_set"].map(lambda p: os.path.basename(p))
res["GO_ID"] = res["Term"].str.extract(r"\((GO:\d+)\)", expand=False)
res.to_csv(_cfg("FIG2_GO_TSV"), sep="\t")
sig = res[res["Adjusted P-value"] <= 0.05]
print(f"terms tested: {len(res):,}; adj. P <= 0.05: {len(sig)} "
      f"({(sig['Gene_set'] == 'human_GO_bp_2025.gmt').sum()} BP)")
print(f"wrote {_cfg('FIG2_NOVEL_GENES_TSV')}\nwrote {_cfg('FIG2_GO_TSV')}")
