#!/usr/bin/env python3
"""
Build a gene-family GFF and metadata table for all paralog families.

Input:  hifiasm_041425_denovoEnhanced_peaks2utr_sorted_perfect.gff
Output: paralog_families.gff       — gene features for all family members (parent + paralogs)
        paralog_families.tsv       — family metadata table
        paralog_families_summary.tsv — per-family summary counts

Paralog suffixes:
  -lN   = Liftoff-transferred paralog (e.g., Cct7-l1)
  -dlN  = de-novo predicted paralog (e.g., Cct7-dl1)
  -rlN  = retrogene-like paralog (e.g., Rbmx-rl1)

For each paralog, the parent name is derived by stripping the -lN/-dlN/-rlN suffix.
"""

import re
import sys
from collections import defaultdict

# Usage: python build_paralog_families.py [input.gff]  (default mirrors config.sh: PARALOG_GFF)
GFF_IN = sys.argv[1] if len(sys.argv) > 1 else "/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj/data/denovo_OctDegus_genome/041425-assembly/hifiasm_041425_denovoEnhanced_peaks2utr_sorted_perfect.gff"
GFF_OUT = "paralog_families.gff"
TSV_OUT = "paralog_families.tsv"
SUMMARY_OUT = "paralog_families_summary.tsv"

# Regex for paralog suffixes
PARALOG_RE = re.compile(r'^(.+)-(dl|rl|l)(\d+)$')

# ---------------------------------------------------------------------------
# Pass 1: collect all gene features (name -> [lines])
# ---------------------------------------------------------------------------
print("Pass 1: Loading gene features from GFF...", file=sys.stderr)

gene_by_name = defaultdict(list)     # Name -> [(attrs, full_gff_line)]
gene_by_id = {}                      # ID   -> (Name, attrs, full_gff_line)

with open(GFF_IN) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue

        attrs = parts[8]
        # Extract ID and Name
        id_m = re.search(r'ID=([^;]+)', attrs)
        name_m = re.search(r'Name=([^;]+)', attrs)

        if not id_m:
            continue

        gid = id_m.group(1)
        name = name_m.group(1) if name_m else gid  # fallback to ID

        gene_by_name[name].append((attrs, line))
        gene_by_id[gid] = (name, attrs, line)

print(f"  Loaded {len(gene_by_id)} gene features ({len(gene_by_name)} unique names)", file=sys.stderr)

# ---------------------------------------------------------------------------
# Pass 2: identify paralogs, derive parent names, build families
# ---------------------------------------------------------------------------
print("Pass 2: Building gene families...", file=sys.stderr)

families = defaultdict(lambda: {"parent_name": None, "parent_entries": [], "paralogs": []})
# families[parent_name] = {parent_name, parent_entries: [(name, gid, paralog_type, copy_num, chr, start, end, strand, length)], paralogs: [...]}

paralog_count = 0
parents_found = 0
parents_missing = 0

for name, entries in gene_by_name.items():
    m = PARALOG_RE.match(name)
    if not m:
        continue  # not a paralog

    parent_name = m.group(1)
    ptype = m.group(2)   # dl, rl, l
    copy_num = int(m.group(3))

    paralog_count += 1

    # Get the first entry's genomic coordinates
    attrs, gff_line = entries[0]
    id_m = re.search(r'ID=([^;]+)', attrs)
    gid = id_m.group(1) if id_m else name
    parts = gff_line.rstrip("\n").split("\t")
    chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
    length = end - start + 1

    paralog_entry = {
        "name": name,
        "gid": gid,
        "ptype": ptype,
        "copy_num": copy_num,
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "length": length,
        "gff_line": gff_line,
    }

    families[parent_name]["paralogs"].append(paralog_entry)

# ---------------------------------------------------------------------------
# Pass 3: find parent genes in the GFF
# ---------------------------------------------------------------------------
print("Pass 3: Matching parent genes...", file=sys.stderr)

for parent_name, fam in families.items():
    if parent_name in gene_by_name:
        entries = gene_by_name[parent_name]
        attrs, gff_line = entries[0]
        id_m = re.search(r'ID=([^;]+)', attrs)
        gid = id_m.group(1) if id_m else parent_name
        parts = gff_line.rstrip("\n").split("\t")
        chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
        length = end - start + 1

        fam["parent_name"] = parent_name
        fam["parent_entries"] = [{
            "name": parent_name,
            "gid": gid,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "length": length,
            "gff_line": gff_line,
        }]
        parents_found += 1
    else:
        fam["parent_name"] = parent_name
        fam["parent_entries"] = []
        parents_missing += 1

print(f"  Families: {len(families)}", file=sys.stderr)
print(f"  Total paralogs: {paralog_count}", file=sys.stderr)
print(f"  Parents found: {parents_found}", file=sys.stderr)
print(f"  Parents missing: {parents_missing}", file=sys.stderr)

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------

# Sort families alphabetically
sorted_families = sorted(families.items(), key=lambda x: x[0])

# --- GFF output ---
print("Writing paralog_families.gff...", file=sys.stderr)
with open(GFF_OUT, "w") as gff_out:
    gff_out.write("##gff-version 3\n")
    gff_out.write(f"##paralog families extracted from {GFF_IN}\n")
    gff_out.write(f"##{len(families)} families, {paralog_count} paralogs\n")
    for parent_name, fam in sorted_families:
        for entry in fam["parent_entries"]:
            gff_out.write(entry["gff_line"])
        for entry in fam["paralogs"]:
            gff_out.write(entry["gff_line"])

# --- TSV metadata ---
print("Writing paralog_families.tsv...", file=sys.stderr)
with open(TSV_OUT, "w") as tsv_out:
    header = [
        "family", "gene_name", "gene_id", "gene_type",
        "chrom", "start", "end", "strand", "length",
        "paralog_type", "copy_num"
    ]
    tsv_out.write("\t".join(header) + "\n")
    for parent_name, fam in sorted_families:
        # Parent first
        for entry in fam["parent_entries"]:
            row = [
                parent_name, entry["name"], entry["gid"], "parent",
                entry["chrom"], str(entry["start"]), str(entry["end"]),
                entry["strand"], str(entry["length"]), "", ""
            ]
            tsv_out.write("\t".join(row) + "\n")
        # Then paralogs
        for entry in fam["paralogs"]:
            row = [
                parent_name, entry["name"], entry["gid"],
                f"paralog_{entry['ptype']}",
                entry["chrom"], str(entry["start"]), str(entry["end"]),
                entry["strand"], str(entry["length"]),
                entry["ptype"], str(entry["copy_num"])
            ]
            tsv_out.write("\t".join(row) + "\n")

# --- Family summary ---
print("Writing paralog_families_summary.tsv...", file=sys.stderr)
with open(SUMMARY_OUT, "w") as sum_out:
    sum_out.write("\t".join([
        "family", "has_parent", "n_paralogs_total",
        "n_paralog_l", "n_paralog_dl", "n_paralog_rl",
        "parent_length", "paralog_total_length"
    ]) + "\n")
    for parent_name, fam in sorted_families:
        has_parent = "yes" if fam["parent_entries"] else "no"
        n_total = len(fam["paralogs"])
        n_l = sum(1 for p in fam["paralogs"] if p["ptype"] == "l")
        n_dl = sum(1 for p in fam["paralogs"] if p["ptype"] == "dl")
        n_rl = sum(1 for p in fam["paralogs"] if p["ptype"] == "rl")
        parent_len = fam["parent_entries"][0]["length"] if fam["parent_entries"] else 0
        paralog_total_len = sum(p["length"] for p in fam["paralogs"])
        sum_out.write("\t".join([
            parent_name, has_parent,
            str(n_total), str(n_l), str(n_dl), str(n_rl),
            str(parent_len), str(paralog_total_len)
        ]) + "\n")

print("Done.", file=sys.stderr)
print(f"  {GFF_OUT}", file=sys.stderr)
print(f"  {TSV_OUT}", file=sys.stderr)
print(f"  {SUMMARY_OUT}", file=sys.stderr)
