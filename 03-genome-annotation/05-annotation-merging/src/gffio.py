"""GFF3 parsing and writing for the annotation merge.

A gene is represented as a :class:`Gene` that keeps the *raw feature rows* it
was parsed from (gene row + its mRNA / exon / CDS / UTR rows), plus the fields
the merge rules need (name, LOC status, CDS coordinates).  Keeping the raw rows
means the merged output preserves the original feature structure of each gene.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Tuple

from .naming import is_loc

# Columns of a GFF3 row (9 fixed columns, attributes in column 9).
SEQID, SOURCE, FTYPE, START, END, SCORE, STRAND, PHASE, ATTRS = range(9)
Row = Tuple[str, str, str, int, int, str, str, str, str]

_ATTR = re.compile(r"^(\w+)=(.+)$")

# Feature types that act as a UTR/CDS parent ("transcript" in the loose sense,
# including Ig/TCR gene segments) and the child features whose Parent points at
# one of them.
TRANSCRIPT_TYPES = {
    "mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA", "snRNA", "snoRNA",
    "C_gene_segment", "V_gene_segment", "D_gene_segment", "J_gene_segment",
}
CHILD_PARENT_TYPES = {"exon", "CDS", "start_codon", "stop_codon",
                      "five_prime_UTR", "three_prime_UTR"}


def parse_attributes(attrs: str) -> dict:
    """Split a GFF3 attribute field into a dict (first value wins per key)."""
    out: dict = {}
    for piece in attrs.split(";"):
        m = _ATTR.match(piece)
        if m:
            key, val = m.group(1), m.group(2)
            out.setdefault(key, val)
    return out


def _first_attr(attrs: str, *keys: str) -> str:
    d = parse_attributes(attrs)
    for k in keys:
        if k in d:
            return d[k]
    return ""


def set_attributes(attrs: str, **updates: str) -> str:
    """Return *attrs* with the given key=value updates applied in place."""
    pieces = attrs.split(";")
    new = []
    done = set()
    for piece in pieces:
        m = _ATTR.match(piece)
        if m and m.group(1) in updates:
            new.append(f"{m.group(1)}={updates[m.group(1)]}")
            done.add(m.group(1))
        else:
            new.append(piece)
    for k, v in updates.items():
        if k not in done:
            new.append(f"{k}={v}")
    return ";".join(x for x in new if x)


@dataclass
class Gene:
    """A gene and all of its child features."""

    seqid: str
    source: str
    start: int
    end: int
    strand: str
    attributes: str
    gene_id: str = ""
    name: str = ""
    biotype: str = ""
    cds_ranges: List[Tuple[int, int]] = field(default_factory=list)
    child_rows: List[Row] = field(default_factory=list)
    # Evidence class used by the naming step:
    #   'replace'  = de-novo gene placed by the replace rule
    #   'rename'   = Liftoff gene renamed by the name-change rule
    #   'de_novo'  = de-novo gene added as a novel gene
    #   'liftoff'  = Liftoff gene kept as-is
    evidence: str = "liftoff"
    replace_concordant: bool = False  # for 'replace': replaced a same-named locus

    @property
    def is_loc(self) -> bool:
        return is_loc(self.name)

    def rows(self) -> List[Row]:
        gene_row = (self.seqid, self.source, "gene", self.start, self.end,
                    ".", self.strand, ".", self.attributes)
        return [gene_row, *self.child_rows]

    def set_name(self, new_name: str) -> None:
        """Set the gene's name, rewriting its ID to ``gene-<new_name>``.

        Used by the name-change rule (Liftoff genes) and when a de-novo gene is
        placed in the output, so that every gene carries an ``ID=gene-*`` that
        downstream tools (and name extraction from ``ID=gene-``) expect.
        """
        self.name = new_name
        self.attributes = set_attributes(self.attributes,
                                         ID=f"gene-{new_name}",
                                         Name=new_name)
        self.gene_id = f"gene-{new_name}"

    def normalize_children(self) -> None:
        """Rebuild ``Parent`` linkage among this gene's children.

        De-novo (Braker) GFF3 files often carry an mRNA whose ``Parent`` still
        points at the old Braker gene ID and CDS/exon rows with no ``Parent`` at
        all.  Re-parent every transcript to the gene and every CDS/exon to the
        transcript whose span contains it; already-valid ``Parent`` values are
        left untouched (the operation is idempotent).
        """
        tx = [(r[START], r[END], _first_attr(r[ATTRS], "ID"))
              for r in self.child_rows
              if r[FTYPE] in TRANSCRIPT_TYPES and _first_attr(r[ATTRS], "ID")]
        if not tx:
            return
        tx_ids = {t for _, _, t in tx}

        def _owner(start: int, end: int) -> Optional[str]:
            for s, e, mid in tx:
                if s <= start and end <= e:
                    return mid
            return tx[0][2] if len(tx) == 1 else None

        rows = []
        for r in self.child_rows:
            ftype, attrs = r[FTYPE], r[ATTRS]
            if ftype in TRANSCRIPT_TYPES and attrs not in ("", "."):
                attrs = set_attributes(attrs, Parent=self.gene_id)
            elif ftype in CHILD_PARENT_TYPES:
                if _first_attr(attrs, "Parent") not in tx_ids:
                    new_par = _owner(r[START], r[END])
                    if new_par:
                        attrs = (f"Parent={new_par}" if attrs in ("", ".")
                                 else set_attributes(attrs, Parent=new_par))
            rows.append((r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], attrs))
        self.child_rows = rows


def parse_gff(path: str) -> List[Gene]:
    """Parse a GFF3/GTF file into a list of :class:`Gene`.

    Features are attached to their owning gene by the ``Parent`` attribute when
    present, otherwise by which gene's span contains the feature.  This is robust
    to files (e.g. ``peaks2utr`` output) where some features precede their gene
    row, which a naive "rows between gene rows" grouping would mis-attribute to
    the *previous* gene.
    """
    _RNA = {"mRNA", "transcript", "lnc_RNA", "ncRNA", "tRNA", "rRNA", "snRNA", "snoRNA"}

    genes: List[Gene] = []
    gene_by_id: dict = {}
    rows: List[Row] = []

    def _norm(gid: str) -> str:
        return gid.removeprefix("gene-").lower()

    # ---- pass 1: collect gene rows ----
    with open(path, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[FTYPE]
            try:
                start, end = int(parts[START]), int(parts[END])
            except ValueError:
                continue
            row = (parts[SEQID], parts[SOURCE], ftype, start, end, parts[SCORE],
                   parts[STRAND], parts[PHASE], parts[ATTRS])
            if ftype == "gene":
                attrs = parts[ATTRS]
                gene_id = _first_attr(attrs, "ID")
                name = _first_attr(attrs, "Name", "gene", "locus_tag")
                if not name and gene_id.startswith("gene-"):
                    name = gene_id[len("gene-"):]
                if not name:
                    name = gene_id
                g = Gene(
                    seqid=parts[SEQID], source=parts[SOURCE],
                    start=start, end=end, strand=parts[STRAND],
                    attributes=attrs, gene_id=gene_id, name=name,
                    biotype=_first_attr(attrs, "gene_biotype", "biotype"),
                )
                genes.append(g)
                gene_by_id[_norm(gene_id)] = g
            else:
                rows.append(row)

    # ---- pass 2: build mRNA -> gene map + a per-seq span index, attach ----
    mrna_to_gene: dict = {}
    for row in rows:
        if row[FTYPE] in _RNA:
            pid = _first_attr(row[ATTRS], "Parent")
            mid = _first_attr(row[ATTRS], "ID")
            if pid:
                mrna_to_gene[mid] = pid

    from collections import defaultdict
    from bisect import bisect_right
    seq_index: dict = defaultdict(list)   # seqid -> [(start, end, gene)]
    for g in genes:
        seq_index[g.seqid].append((g.start, g.end, g))
    seq_starts: dict = {}
    for seq, arr in seq_index.items():
        arr.sort(key=lambda x: x[0])
        seq_starts[seq] = [x[0] for x in arr]

    for row in rows:
        g = _gene_for(row, gene_by_id, mrna_to_gene, seq_index, seq_starts)
        if g is None:
            continue
        g.child_rows.append(row)
        if row[FTYPE] == "CDS":
            g.cds_ranges.append((row[START], row[END]))
        if row[START] < g.start:
            g.start = row[START]
        if row[END] > g.end:
            g.end = row[END]

    return genes


def _gene_for(row: Row, gene_by_id: dict, mrna_to_gene: dict,
              seq_index: dict, seq_starts: dict) -> Optional[Gene]:
    """Resolve the gene that owns a feature row."""
    pid = _first_attr(row[ATTRS], "Parent")
    if pid:
        gid = pid
        if gid.startswith("rna-") or (gid.lower() in mrna_to_gene):
            gid = mrna_to_gene.get(pid, mrna_to_gene.get(gid, ""))
        g = gene_by_id.get(gid.removeprefix("gene-").lower())
        if g is not None:
            return g
    # fallback: the gene whose span contains the feature start (binary search)
    from bisect import bisect_right
    arr = seq_index.get(row[SEQID])
    if not arr:
        return None
    i = bisect_right(seq_starts[row[SEQID]], row[START]) - 1
    while i >= 0:
        s, e, g = arr[i]
        if e < row[START]:
            break
        if s <= row[START]:
            return g
        i -= 1
    return None


def write_gff(genes: Sequence[Gene], path: str, header: Iterable[str] = ()) -> None:
    """Write *genes* (sorted by seqid, start) to *path* as GFF3."""
    ordered = sorted(genes, key=lambda g: (g.seqid, g.start, g.end, g.strand))
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for h in header:
            fh.write(f"##{h}\n" if not h.startswith("#") else f"{h}\n")
        for gene in ordered:
            for row in gene.rows():
                fh.write("\t".join(str(x) for x in row) + "\n")
