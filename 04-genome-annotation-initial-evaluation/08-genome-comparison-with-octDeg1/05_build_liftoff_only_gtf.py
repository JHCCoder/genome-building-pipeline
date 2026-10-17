#!/usr/bin/env python3
"""
Build a Cell Ranger ARC GTF for OctDeg2.0 that carries ONLY the OctDeg1.0 RefSeq gene
models lifted over by Liftoff (no BRAKER genes, no peaks2utr UTRs), restricted to the
same transcripts as the OctDeg1.0 Cell Ranger reference. Panel D then compares the two
assemblies with matched gene models, so a difference reflects contiguity, not annotation.

Matching:
  * transcripts with an NCBI accession (XM_/XR_/NM_) are kept when the same accession is
    in the OctDeg1.0 reference GTF;
  * reference transcripts without an accession (Cell Ranger's "unassigned_transcript_N",
    mostly tRNAs) are matched by gene: a Liftoff transcript is kept when its gene symbol
    is one of those reference genes.
  * gene_id / gene_name / biotypes are taken from the OctDeg1.0 reference, so both
    references use identical gene identifiers.
Reference transcripts that Liftoff did not place on OctDeg2.0 are simply absent; their
count is printed and written to the summary file.

Output records are grouped gene -> transcript -> exon/CDS, which is the order
cellranger-arc mkref requires (each exon must follow its transcript row).

Inputs and outputs come from config.sh:
  FIG3D_OCTDEG1_CR_REF   Cell Ranger ARC OctDeg1.0 reference (genes/genes.gtf.gz)
  FIG3D_LIFTOFF_GFF      Liftoff GFF of OctDeg1.0 RefSeq onto MITO_ASSEMBLY
  FIG3D_LIFTOFF_GTF      output GTF (a *.summary.txt is written next to it)

Run with the python-visualizations env:
    python 05_build_liftoff_only_gtf.py
"""

import gzip
import os
import re
import subprocess
import sys
from collections import Counter, OrderedDict, defaultdict


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


REF_GTF = os.path.join(_cfg("FIG3D_OCTDEG1_CR_REF"), "genes", "genes.gtf.gz")
LIFTOFF_GFF = _cfg("FIG3D_LIFTOFF_GFF")
OUT_GTF = _cfg("FIG3D_LIFTOFF_GTF")
if os.path.exists(OUT_GTF):
    sys.exit(f"ERROR: {OUT_GTF} exists; move it aside (e.g. *.pre-<date>) before rebuilding")
os.makedirs(os.path.dirname(OUT_GTF), exist_ok=True)

# --- 1. OctDeg1.0 reference transcripts -------------------------------------------------
GTF_ATTR = re.compile(r'(\S+) "([^"]*)"')
ref_tx = {}                       # transcript_id -> attrs
ref_genes = set()
for line in gzip.open(REF_GTF, "rt"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[2] != "transcript":
        continue
    a = dict(GTF_ATTR.findall(f[8]))
    ref_tx[a["transcript_id"]] = a
    ref_genes.add(a["gene_id"])
unassigned_genes = {a["gene_id"] for t, a in ref_tx.items() if t.startswith("unassigned_transcript_")}
ref_gene_attrs = {}
for a in ref_tx.values():
    ref_gene_attrs.setdefault(a["gene_id"], a)

# --- 2. Liftoff GFF: collect transcripts and their exon/CDS children --------------------
def gff_attrs(s):
    return dict(kv.split("=", 1) for kv in s.split(";") if "=" in kv)

tx_row = OrderedDict()            # liftoff transcript ID -> GFF fields
children = defaultdict(list)      # liftoff transcript ID -> [GFF fields of exon/CDS]
for line in open(LIFTOFF_GFF):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    a = gff_attrs(f[8])
    par = a.get("Parent", "")
    if f[2] in ("exon", "CDS"):
        for p in par.split(","):
            children[p].append(f)
    elif par.startswith("gene-") and "ID" in a:
        tx_row[a["ID"]] = f

# --- 3. Select transcripts matching the reference --------------------------------------
kept = []                         # (liftoff ID, transcript_id, ref attrs)
for lid, f in tx_row.items():
    a = gff_attrs(f[8])
    acc = lid[4:] if lid.startswith("rna-") else lid
    sym = a["Parent"][5:]
    if acc in ref_tx:
        kept.append((lid, acc, ref_tx[acc]))
    elif sym in unassigned_genes:
        kept.append((lid, acc, ref_gene_attrs[sym]))

# --- 4. Write GTF grouped gene -> transcript -> exon/CDS --------------------------------
def gtf_attr(d):
    return " ".join(f'{k} "{v}";' for k, v in d.items())

by_gene = OrderedDict()
for lid, tid, ra in kept:
    by_gene.setdefault(ra["gene_id"], []).append((lid, tid, ra))

n_exon = n_cds = 0
no_exon = []
with open(OUT_GTF, "w") as out:
    out.write("#OctDeg2.0 Liftoff-only GTF (OctDeg1.0 RefSeq models), matched to the OctDeg1.0 "
              "Cell Ranger reference; built by 05_build_liftoff_only_gtf.py\n")
    for gid, txs in by_gene.items():
        rows = [tx_row[lid] for lid, _, _ in txs]
        chrom, strand = rows[0][0], rows[0][6]
        start, end = min(int(r[3]) for r in rows), max(int(r[4]) for r in rows)
        ra0 = txs[0][2]
        gattr = OrderedDict(gene_id=gid, gene_name=ra0.get("gene", gid))
        if ra0.get("gene_biotype"):
            gattr["gene_biotype"] = ra0["gene_biotype"]
        out.write("\t".join([chrom, "Liftoff", "gene", str(start), str(end), ".", strand, ".",
                             gtf_attr(gattr)]) + "\n")
        for lid, tid, ra in txs:
            r = tx_row[lid]
            exons = sorted((c for c in children[lid] if c[2] == "exon"), key=lambda c: int(c[3]))
            if not exons:
                no_exon.append(tid)
                continue
            tattr = OrderedDict(gattr)
            tattr["transcript_id"] = tid
            if ra.get("transcript_biotype"):
                tattr["transcript_biotype"] = ra["transcript_biotype"]
            out.write("\t".join([r[0], "Liftoff", "transcript", r[3], r[4], ".", r[6], ".",
                                 gtf_attr(tattr)]) + "\n")
            for c in exons:
                out.write("\t".join([c[0], "Liftoff", "exon", c[3], c[4], ".", c[6], ".",
                                     gtf_attr(tattr)]) + "\n")
                n_exon += 1
            for c in sorted((c for c in children[lid] if c[2] == "CDS"), key=lambda c: int(c[3])):
                out.write("\t".join([c[0], "Liftoff", "CDS", c[3], c[4], ".", c[6], c[7],
                                     gtf_attr(tattr)]) + "\n")
                n_cds += 1

kept_tx = {tid for _, tid, _ in kept} - set(no_exon)
kept_genes = {ra["gene_id"] for _, tid, ra in kept if tid in kept_tx}
missing_genes = sorted(ref_genes - kept_genes)
biotype_missing = Counter(ref_gene_attrs[g].get("transcript_biotype", "NA") for g in missing_genes)
summary = [
    f"OctDeg1.0 reference : {len(ref_tx):,} transcripts, {len(ref_genes):,} genes ({REF_GTF})",
    f"Liftoff GFF         : {len(tx_row):,} transcripts ({LIFTOFF_GFF})",
    f"Kept for OctDeg2.0  : {len(kept_tx):,} transcripts, {len(kept_genes):,} genes, "
    f"{n_exon:,} exons, {n_cds:,} CDS",
    f"Reference genes not placed by Liftoff: {len(missing_genes):,} "
    f"({100 * len(missing_genes) / len(ref_genes):.1f}%); by biotype: {dict(biotype_missing.most_common())}",
    f"Kept transcripts dropped for lacking exons: {len(no_exon):,}",
    f"Output: {OUT_GTF}",
]
print("\n".join(summary))
with open(OUT_GTF.rsplit(".", 1)[0] + ".summary.txt", "w") as fh:
    fh.write("\n".join(summary) + "\n\nReference genes not placed by Liftoff:\n" + "\n".join(missing_genes) + "\n")
