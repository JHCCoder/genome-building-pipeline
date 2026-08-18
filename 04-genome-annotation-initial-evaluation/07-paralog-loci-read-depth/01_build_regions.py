#!/usr/bin/env python3
"""
Build region BED files for the paralog-vs-single-copy HiFi read-coverage
analysis (Supplemental Figure S7).

Classifies annotated genes into four tiers plus a BUSCO control set:

  T1 (G2): paralog COPIES (name ends in -lN/-dlN/-rlN), NO segdup overlap
  T2 (G4): paralog COPIES overlapping segmental duplications (segdups)
  T3 (G3): non-paralog genes overlapping segdups
  T4 (G1): non-paralog genes NOT overlapping segdups
  BUSCO:   single-copy orthologs (glires_odb10)

"Paralog" follows the strict copy definition of build_paralog_families.py
(04-…/03-paralog-assessment-by-read-depth): a gene is a paralog only if its
name carries a -lN / -dlN / -rlN copy suffix (the parent gene is NOT itself
a paralog copy).  All chromosomes are kept (no chr1-28 restriction).

Each locus region = gene body +/- FLANK_BP (default 5 kb).

Outputs (written to this script's directory):
  regions.tsv               - full gene-level tier classification table
  regions.bed               - gene +/- flank regions (chr, start, end, gene, tier, overlaps_dup)
  gene_body.bed             - gene-body-only regions (no flanks)
  busco_regions.bed         - BUSCO gene bodies
  busco_regions_5kb_flank.bed - BUSCO gene bodies +/- flank

Run with a conda env that has pandas/numpy (e.g. `genome-assembly`):

    conda activate genome-assembly
    python 01_build_regions.py
"""

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# RUN PARAMETERS
# ---------------------------------------------------------------------------
FLANK_BP = 5000            # bp of upstream/downstream flank added to each gene

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
# Project root — edit for your environment.
PROJ_ROOT = "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"

# This script's directory (all outputs are written here).
HERE = Path(__file__).resolve().parent

# Gene x segdup overlap table, produced from the NEW merged annotation
# (hifiasm-041425-denovoEnhanced_peaks2utr_sorted.agat.gff3) by
# 06-segmental-duplication-analysis/04_segdup_gene_overlap.ipynb.
# Re-run that notebook after a new annotation before running this script.
GENE_FILE = (
    f"{PROJ_ROOT}/code/github-code-to-share/06-segmental-duplication-analysis/"
    "hifiasm_gene_segDup_overlapInfo_081726.tsv"
)

# Chromosome lengths for the final assembly (one "<chrom>\t<len>" pair per line).
GENOME_LENGTH_FILE = (
    f"{PROJ_ROOT}/code/command-line-script/contig-coverage/"
    "hifiasm_041425_scaffolded_juiceBox_sorted_hardMasked_chrAssigned_genome_length.txt"
)

# BUSCO run full_table.tsv (glires_odb10); assembly-level, annotation-independent.
BUSCO_FILE = (
    f"{PROJ_ROOT}/output/outputs-from-busco-ortholog-alignment/"
    "hifiasm-041425-primary-scaffolded/run_glires_odb10/full_table.tsv"
)

OUT_DIR = str(HERE)

# ============================================================
# Load data
# ============================================================
print("Loading gene-segdup overlap file...")
genes = pd.read_csv(GENE_FILE, sep="\t")
print(f"  Total genes loaded: {len(genes)}")

print("Loading chromosome lengths...")
chr_lengths = {}
with open(GENOME_LENGTH_FILE) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            chr_lengths[parts[0]] = int(parts[1])
print(f"  Loaded {len(chr_lengths)} chromosome lengths")

# ============================================================
# Check chromosomes against the length file (all chromosomes kept)
# ============================================================
missing_chroms = set(genes["Chromosome"].unique()) - set(chr_lengths.keys())
if missing_chroms:
    print(f"  WARNING: Chromosomes missing from length file: {missing_chroms}")

# ============================================================
# Extract base gene name and classify paralogs
# ============================================================
print("Extracting base gene names and classifying paralogs...")


# Strict paralog-copy suffix — same rule as build_paralog_families.py in
# 04-…/03-paralog-assessment-by-read-depth: -lN / -dlN / -rlN only.
# The overlap file's `gene` column is UPPERCASED, so match case-insensitively.
PARALOG_COPY_RE = re.compile(r"^(.+)-(dl|rl|l)\d+$", re.IGNORECASE)


def get_parent_name(name):
    """Return the parent gene name (strip a trailing -lN/-dlN/-rlN suffix)."""
    if pd.isna(name):
        return name
    m = PARALOG_COPY_RE.match(str(name))
    return m.group(1) if m else str(name)


def is_paralog_copy(name):
    """True if the gene is a paralog COPY (name ends in -lN/-dlN/-rlN)."""
    return bool(PARALOG_COPY_RE.match(str(name))) if pd.notna(name) else False


genes["base_gene"] = genes["gene"].apply(get_parent_name)
genes["is_paralog"] = genes["gene"].apply(is_paralog_copy)

# Family size = number of genes sharing a parent name (parent + its copies).
base_counts = genes["base_gene"].value_counts()
genes["num_copies"] = genes["base_gene"].map(base_counts)

# Families = unique parent names among the paralog copies (matches
# build_paralog_families.py's family count).
paralog_bases = set(genes.loc[genes["is_paralog"], "base_gene"])

# ============================================================
# Classify into tiers
# ============================================================
print("Classifying into tiers...")


def classify_tier(row):
    if row["is_paralog"]:
        return "T2" if row["overlaps_dup"] else "T1"
    return "T3" if row["overlaps_dup"] else "T4"


genes["tier"] = genes.apply(classify_tier, axis=1)

# ============================================================
# Compute region coordinates (+/- flank)
# ============================================================
print(f"Computing regions with +/-{FLANK_BP} bp flanks...")


def compute_region(row):
    """Region start/end with flanks, clamped to chromosome boundaries."""
    chrom = row["Chromosome"]
    chr_len = chr_lengths.get(chrom, None)

    gene_start = int(row["Start"])  # 0-based
    gene_end = int(row["End"])      # 1-based exclusive

    region_start = max(0, gene_start - FLANK_BP)
    region_end = min(chr_len, gene_end + FLANK_BP) if chr_len is not None else gene_end + FLANK_BP

    return pd.Series([region_start, region_end, gene_start, gene_end])


genes[["region_start", "region_end", "gene_start_bed", "gene_end_bed"]] = genes.apply(
    compute_region, axis=1
)

genes["region_length"] = genes["region_end"] - genes["region_start"]
genes["gene_length"] = genes["gene_end_bed"] - genes["gene_start_bed"]

# ============================================================
# Summary statistics
# ============================================================
print("\n" + "=" * 60)
print("TIER SUMMARY")
print("=" * 60)
for tier in ["T1", "T2", "T3", "T4"]:
    subset = genes[genes["tier"] == tier]
    print(f"  {tier}: {len(subset):>6} genes  "
          f"(median gene len: {subset['gene_length'].median():,.0f} bp, "
          f"median region len: {subset['region_length'].median():,.0f} bp)")

print(f"\n  Total: {len(genes)} genes")
print(f"  Paralog families: {len(paralog_bases)}")
print(f"  Paralog copies (non-parent): {int(genes['is_paralog'].sum())}")
print(f"  Non-paralog genes (singletons + parents): {int((~genes['is_paralog']).sum())}")
print(f"  Overlap segdup: {int(genes['overlaps_dup'].sum())}")

# ============================================================
# Write outputs
# ============================================================
print("\nWriting outputs...")

tsv_out = os.path.join(OUT_DIR, "regions.tsv")
cols_out = [
    "gene", "base_gene", "Chromosome", "gene_start_bed", "gene_end_bed",
    "region_start", "region_end", "region_length", "gene_length",
    "tier", "is_paralog", "num_copies", "overlaps_dup",
    "Strand", "source",
]
genes[cols_out].to_csv(tsv_out, sep="\t", index=False)
print(f"  {tsv_out}")

bed_out = os.path.join(OUT_DIR, "regions.bed")
bed_extended = pd.DataFrame({
    "chrom": genes["Chromosome"],
    "start": genes["region_start"].astype(int),
    "end": genes["region_end"].astype(int),
    "gene": genes["gene"],
    "tier": genes["tier"],
    "overlaps_dup": genes["overlaps_dup"].astype(str),
})
bed_extended.to_csv(bed_out, sep="\t", index=False, header=False)
print(f"  {bed_out} (cols: chrom, start, end, gene, tier, overlaps_dup)")

bed_gene_out = os.path.join(OUT_DIR, "gene_body.bed")
bed_gene = pd.DataFrame({
    "chrom": genes["Chromosome"],
    "start": genes["gene_start_bed"].astype(int),
    "end": genes["gene_end_bed"].astype(int),
    "gene": genes["gene"],
    "tier": genes["tier"],
    "overlaps_dup": genes["overlaps_dup"].astype(str),
})
bed_gene.to_csv(bed_gene_out, sep="\t", index=False, header=False)
print(f"  {bed_gene_out} (gene body only, no flanks)")

# BUSCO single-copy genes as an additional control set
print("\nLoading BUSCO single-copy genes...")
# BUSCO full_table.tsv has 2 comment lines, then a header that also starts with "# ".
# Skip the first 2 lines, use line 3 as header, strip "# " from column names.
busco = pd.read_csv(BUSCO_FILE, sep="\t", skiprows=2, header=0)
busco.columns = [c.lstrip("# ") for c in busco.columns]
busco_complete = busco[busco["Status"] == "Complete"].copy()
print(f"  BUSCO Complete (single-copy): {len(busco_complete)}")

# BUSCO uses 1-based coordinates; convert to 0-based BED. For minus-strand
# genes, Gene Start > Gene End, so take min/max.
busco_gs = busco_complete["Gene Start"].astype(int)
busco_ge = busco_complete["Gene End"].astype(int)
busco_bed_start = np.minimum(busco_gs, busco_ge) - 1  # 1-based -> 0-based
busco_bed_end = np.maximum(busco_gs, busco_ge)

busco_bed = pd.DataFrame({
    "chrom": busco_complete["Sequence"],
    "start": busco_bed_start,
    "end": busco_bed_end,
    "gene": busco_complete["Busco id"],
    "tier": "BUSCO",
    "overlaps_dup": "False",
})
busco_bed_out = os.path.join(OUT_DIR, "busco_regions.bed")
busco_bed.to_csv(busco_bed_out, sep="\t", index=False, header=False)
print(f"  {busco_bed_out} ({len(busco_bed)} BUSCO regions)")

busco_flank_bed = pd.DataFrame({
    "chrom": busco_complete["Sequence"],
    "start": (busco_bed_start - FLANK_BP).clip(lower=0),
    "end": busco_bed_end + FLANK_BP,
    "gene": busco_complete["Busco id"],
    "tier": "BUSCO",
    "overlaps_dup": "False",
})
busco_flank_out = os.path.join(OUT_DIR, "busco_regions_5kb_flank.bed")
busco_flank_bed.to_csv(busco_flank_out, sep="\t", index=False, header=False)
print(f"  {busco_flank_out} ({len(busco_flank_bed)} BUSCO regions +/-{FLANK_BP}bp)")

print("\nDone!")
