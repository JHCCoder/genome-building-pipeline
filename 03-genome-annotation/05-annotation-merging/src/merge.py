"""Core merge logic: the three enhancement rules.

Rules (see README.md for rationale):

1. **name-change** — optionally rename Liftoff genes using a user-provided
   mapping (e.g. correcting an unannotated ``LOC`` name from a comparative
   analysis).  Liftoff structure is kept; only the name changes.

2. **replace** — a functionally-named de-novo gene replaces a Liftoff gene it
   overlaps **only** when the Liftoff gene is an unannotated ``LOC`` or when the
   de-novo gene *spans multiple transferred genes* and one of them shares its
   name (the de-novo gene is re-annotating a locus that Liftoff split apart):
     * Liftoff gene is ``LOC##########`` → replaced (LOC rule);
     * de-novo name == a Liftoff name (case/suffix-insensitive) **and** the
       de-novo gene overlaps ≥ 2 Liftoff genes → the concordant Liftoff gene is
       replaced (concordant rule) — the de-novo gene is a better, unified model
       of a fragmented locus;
     * de-novo name == a single overlapping Liftoff gene (a 1:1 concordant
       overlap) → the Liftoff gene is **kept** (concordant-1:1 rule): they are
       the same gene and there is no reason to drop the reference structure;
     * otherwise (named Liftoff gene with a different name, e.g. a paralog
       locus) → the Liftoff gene is **kept** (discordant rule), because the
       reference-derived annotation should not be silently discarded for a
       de-novo prediction carrying a different name.
   The decision is made **per overlapping Liftoff gene**, so a single de-novo
   gene can replace the concordant member of a tandem pair while the
   differently-named neighbour is retained.

3. **add** — a functionally-named de-novo gene is added as a novel gene when its
   CDS does not overlap any retained Liftoff gene's CDS (strict by default),
   it does not duplicate an existing gene name, and it passes optional
   expression / CDS-length filters.
"""
from __future__ import annotations

import re

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .gffio import Gene, set_attributes
from .naming import (base_name, base_name_aliases, has_suffix, names_concordant,
                     normalize, suffix_type)
from .overlap import IntervalIndex, cds_overlap_fraction, overlap_fraction


# RefSeq paralog locus tags: a gene ID like 'gene-Iqcn-2' marks a reference-
# identified paralog copy whose true name is the ID-derived locus (Iqcn-2),
# distinct from the shared gene symbol (Iqcn).  Such names are preserved.
_REF_PARALOG_ID = re.compile(r"^gene-.+-\d+$", re.IGNORECASE)


def _ref_paralog_name(g: Gene) -> str:
    """Return the reference locus name (e.g. 'Iqcn-2') if *g* carries a RefSeq
    paralog locus tag, else ''."""
    gid = g.gene_id or ""
    if _REF_PARALOG_ID.match(gid):
        return gid[len("gene-"):]
    return ""


@dataclass
class MergeConfig:
    require_name: bool = True                 # skip de-novo genes with no functional name
    min_cds_bp: int = 0                       # min de-novo CDS length to consider
    overlap_fraction: float = 0.0             # min gene-span overlap (fraction) to trigger replace
    cds_overlap_fraction: float = 0.0         # max CDS overlap (fraction) allowed for a novel add
    skip_duplicate_base_name: bool = True     # don't add a de-novo gene whose base name already exists
    rename_map: Dict[str, str] = field(default_factory=dict)  # name/gene_id -> new name
    bigwig: Optional[str] = None              # optional expression track (BigWig) for the add rule
    min_expression: float = 0.0               # min mean signal to keep a de-novo gene in the add rule
    bigwig_chrom_style: str = "ucsc"          # 'ucsc' (chr1) or 'ensembl' (1)


@dataclass
class Decision:
    """A per-gene record of what the merge did."""
    category: str
    liftoff_gene: str
    denovo_gene: str
    detail: str = ""


@dataclass
class MergeResult:
    final_genes: List[Gene]
    decisions: List[Decision]

    @property
    def counts(self) -> Dict[str, int]:
        c: Dict[str, int] = {}
        for d in self.decisions:
            c[d.category] = c.get(d.category, 0) + 1
        return c


def _cds_length(g: Gene) -> int:
    return sum(e - s + 1 for s, e in g.cds_ranges)


# Gene symbols conventionally written all-caps (as in the reference / original
# annotation): mitochondrial genes, complement components and coagulation factors.
_ALL_CAPS_SYMBOLS = {
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB",
    "C2", "C3", "C5", "C6", "C7", "C9",
    "F2", "F3", "F5", "F7", "F8", "F9", "F10", "F11", "F12",
}


def _canonical_base(base: str) -> str:
    """Match the original annotation's symbol case: all-caps for the
    conventional set and LOC genes, title case otherwise."""
    if re.match(r"^LOC\d{9}$", base, re.IGNORECASE):
        return base.upper()
    up = base.upper()
    if up in _ALL_CAPS_SYMBOLS:
        return up
    return base[:1].upper() + base[1:].lower()


def _canonical_name(name: str) -> str:
    """Match the original annotation's gene-name case: title-case symbol with a
    lowercase paralog suffix (-l1 / -dl1 / -rl1)."""
    m = re.match(r"^(.+)-([dD]?[rR]?[lL]\d+)$", name)
    if m:
        return _canonical_base(m.group(1)) + "-" + m.group(2).lower()
    return _canonical_base(name)


def assign_names(genes: List[Gene]) -> List[Gene]:
    """Assign unique gene names following the evidence-ranked convention.

    Within a gene family (same suffix-stripped base name), the gene that keeps
    the **plain parent name** is chosen by evidence priority:

      1. a **replace**-introduced de-novo gene whose replacement was *concordant*
         (it took over the reference locus that carried the same name);
      2. a gene **renamed** by the name-change rule that holds the plain name;
      3. a kept **Liftoff** copy with the plain name;
      4. a **de-novo** (added) gene.

    Every other copy is a paralog and is suffixed by its evidence class —
    ``-RL`` for replace copies, ``-L`` for renamed paralogs, ``-DL`` for
    de-novo copies, and an ``-L`` suffix for kept Liftoff paralogs — numbered
    by genomic order.  An already-correct suffix is preserved when it is free.
    Returns the (mutated) gene list.
    """
    from collections import defaultdict

    groups: Dict[str, List[Gene]] = defaultdict(list)
    for g in genes:
        groups[base_name(g.name)].append(g)

    used: set = set()

    def _assign(g: Gene, name: str) -> None:
        g.set_name(_canonical_name(name))   # match the original annotation's case
        used.add(name.upper())

    def _suffix(g: Gene, base: str, suffix: str) -> None:
        # number sequentially (copies are processed in genomic order per type),
        # so a lone paralog is always -L1 / -DL1 / -RL1 — no gaps.
        n = 1
        while True:
            cand = f"{base}-{suffix}{n}"
            if cand.upper() not in used:
                _assign(g, cand)
                return
            n += 1

    def _plain(g: Gene) -> bool:
        return not has_suffix(g.name)

    def _stype(g: Gene) -> str:
        ev = g.evidence or "liftoff"
        if ev == "replace":
            return "RL"
        if ev == "rename":
            return "L"
        if ev == "de_novo":
            return "DL"
        return suffix_type(g.name) or "L"

    for members in groups.values():
        members.sort(key=lambda g: (g.seqid, g.start))
        # preserve RefSeq paralog locus names (e.g. Iqcn-2) when they are unique;
        # reference genes with genuinely duplicate names (e.g. two H1-0 loci) are
        # still disambiguated below.
        ref_names: Dict[int, str] = {}
        for g in members:
            rn = _ref_paralog_name(g)
            if rn and rn.upper() not in used:
                _assign(g, rn)
                ref_names[id(g)] = rn
        renamable = [g for g in members if id(g) not in ref_names]
        if not renamable:
            continue
        parent = (next((g for g in renamable if g.evidence == "replace" and g.replace_concordant), None)
                  or next((g for g in renamable if g.evidence == "rename" and _plain(g)), None)
                  or next((g for g in renamable if g.evidence == "liftoff" and _plain(g)), None)
                  or next((g for g in renamable if g.evidence == "de_novo" and _plain(g)), None)
                  or renamable[0])
        base = base_name(parent.name)
        if not has_suffix(parent.name) and parent.name.upper() not in used:
            _assign(parent, parent.name)   # keep the original (reference) name/case if free
        elif base.upper() not in used:
            _assign(parent, base)          # parent keeps the plain family name
        else:
            _suffix(parent, base, _stype(parent))  # name taken -> re-suffix like a copy
        by_type: Dict[str, List[Gene]] = defaultdict(list)
        for g in renamable:
            if g is parent:
                continue
            by_type[_stype(g)].append(g)
        for t, copies in by_type.items():  # renamable already sorted by (seqid, start)
            for g in copies:
                _suffix(g, base, t)
    # record each gene's origin so downstream figures can categorize by it
    for g in genes:
        g.attributes = set_attributes(g.attributes, Origin=g.evidence or "liftoff")
    return genes


def _expression_mean(g: Gene, cfg: MergeConfig) -> float:
    """Mean BigWig signal over a gene's CDS (requires pyBigWig)."""
    import pyBigWig  # imported lazily; only needed when cfg.bigwig is set

    bw = pyBigWig.open(cfg.bigwig)
    try:
        chrom = g.seqid
        if cfg.bigwig_chrom_style == "ensembl" and not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        if cfg.bigwig_chrom_style == "ucsc" and chrom.startswith("chr"):
            chrom = chrom
        total = 0.0
        bases = 0
        for s, e in g.cds_ranges:
            try:
                vals = bw.values(chrom, s - 1, e)
            except RuntimeError:
                continue
            if vals:
                total += sum(v for v in vals if v is not None and v == v)
                bases += len([v for v in vals if v is not None and v == v])
        return total / bases if bases else 0.0
    finally:
        bw.close()


def _gene_identifiers(g: Gene) -> set:
    """All identifiers a Liftoff gene can be referred to by."""
    ids = {normalize(g.name), normalize(g.gene_id),
           normalize(g.gene_id).lstrip("GENE-")}
    from .gffio import parse_attributes
    for part in parse_attributes(g.attributes).get("Dbxref", "").split(","):
        if part.startswith("GeneID:"):
            ids.add(part[len("GeneID:"):].strip().upper())
    return ids


def _find_liftoff(liftoff_genes: List[Gene], ident: str) -> Optional[Gene]:
    """Find a Liftoff gene by name, gene-ID (with/without 'gene-') or GeneID.

    ``gene-LOC111813055`` (a GeneID-derived identifier) matches the gene whose
    ``Dbxref=GeneID:111813055`` even if its ``Name=`` is ``Rnf103``.
    """
    key = normalize(ident).lstrip("GENE-")
    geneid_key = key[3:] if key.startswith("LOC") else key
    for g in liftoff_genes:
        ids = _gene_identifiers(g)
        if key in ids or geneid_key in ids:
            return g
    return None


def _find_denovo(denovo_genes: List[Gene], ident: str,
                 near: Optional[Gene]) -> Optional[Gene]:
    """Find a de-novo gene by gene-ID or by base-name (preferring one that
    overlaps *near*).  Multi-name ``Name=`` attributes are matched on any name."""
    gid = normalize(ident)
    exact = [d for d in denovo_genes if normalize(d.gene_id) == gid]
    if exact:
        return exact[0]
    target = base_name_aliases(ident)
    by_name = [d for d in denovo_genes if base_name_aliases(d.name) & target]
    if not by_name:
        return None
    if near is not None:
        for d in by_name:
            if (d.seqid == near.seqid and d.start <= near.end and d.end >= near.start):
                return d
    return by_name[0]


def merge_with_lists(liftoff_genes: List[Gene],
                     denovo_genes: List[Gene],
                     cfg: MergeConfig,
                     rename_map: Optional[Dict[str, str]] = None,
                     add_list: Optional[Dict[str, str]] = None,
                     replace_list: Optional[List[Tuple[str, str]]] = None,
                     ) -> MergeResult:
    """List-driven merge reproducing the original curated logic.

    * ``rename_map``    {liftoff ident: new name}  -> name-change rule
    * ``add_list``      {de-novo gene id: new name} -> add rule
    * ``replace_list``  [(liftoff ident, new name)] -> replace rule

    The replace rule applies the **fix** that protects reference-derived genes:
    a named (non-LOC) Liftoff gene is only replaced when the new name is
    concordant with its own name; otherwise it is kept.
    """
    decisions: List[Decision] = []
    rename_map = rename_map or {}
    add_list = add_list or {}
    replace_list = replace_list or []

    # ---- name-change rule ----------------------------------------------------
    for g in liftoff_genes:
        key = normalize(g.name) or normalize(g.gene_id).lstrip("GENE-")
        if key in rename_map:
            new_name = rename_map[key]
            g.set_name(new_name)
            g.evidence = "rename"
            decisions.append(Decision("renamed", g.gene_id or g.name, "", f"-> {new_name}"))

    removed: set = set()
    placed: set = set()

    # ---- replace rule (with the discordant-named fix) -------------------------
    for lo_ident, new_name in replace_list:
        l = _find_liftoff(liftoff_genes, lo_ident)
        if l is None or id(l) in removed:
            continue
        if not l.is_loc and not names_concordant(l.name, new_name):
            # FIX: keep reference-derived named genes when the replacing name differs
            decisions.append(Decision(
                "replace_skipped_discordant_kept", l.name, new_name,
                f"{l.seqid}:{l.start}-{l.end} kept (name differs)"))
            continue
        d = _find_denovo(denovo_genes, new_name, near=l)
        if d is None:
            decisions.append(Decision("replace_skipped_no_denovo", l.name, new_name, ""))
            continue
        removed.add(id(l))
        placed.add(id(d))
        d.set_name(new_name)   # apply the list's intended output name (e.g. CLEC4A-rl1)
        d.evidence = "replace"
        if not l.is_loc and names_concordant(l.name, new_name):
            d.replace_concordant = True   # replaced a same-named locus -> parent candidate
        category = "replaced_loc" if l.is_loc else "replaced_concordant"
        decisions.append(Decision(
            category, l.name, new_name,
            f"{l.seqid}:{l.start}-{l.end} <- de-novo {d.gene_id}"))

    # ---- add rule ------------------------------------------------------------
    for gid, new_name in add_list.items():
        d = _find_denovo(denovo_genes, gid, near=None)
        if d is None or id(d) in placed:
            continue
        d.set_name(new_name or d.name)
        d.evidence = "de_novo"
        placed.add(id(d))
        decisions.append(Decision(
            "added", "", d.name, f"{d.seqid}:{d.start}-{d.end} (de-novo {d.gene_id})"))

    # ---- assemble -------------------------------------------------------------
    final_genes: List[Gene] = [g for g in liftoff_genes if id(g) not in removed]
    for d in denovo_genes:
        if id(d) in placed:
            d.set_name(d.name)          # ID=gene-<name>;Name=<name> for every placed gene
            final_genes.append(d)
    assign_names(final_genes)
    for g in final_genes:
        g.normalize_children()  # rebuild Parent linkage (de-novo genes carry stale/missing Parent)
    return MergeResult(final_genes=final_genes, decisions=decisions)


def merge_annotations(liftoff_genes: List[Gene],
                      denovo_genes: List[Gene],
                      cfg: MergeConfig) -> MergeResult:
    decisions: List[Decision] = []

    # ---- rule 1: name-change (optional mapping) ----------------------------
    for g in liftoff_genes:
        key = normalize(g.name) or normalize(g.gene_id).lstrip("GENE-")
        if key in cfg.rename_map:
            new_name = cfg.rename_map[key]
            if new_name and not names_concordant(g.name, new_name):
                g.set_name(new_name)
                g.evidence = "rename"
                decisions.append(Decision("renamed", g.gene_id or g.name, "", f"-> {new_name}"))

    # ---- build candidate list from the de-novo annotation -------------------
    candidates: List[Gene] = []
    for d in denovo_genes:
        if cfg.require_name and not normalize(d.name):
            decisions.append(Decision("denovo_skipped_no_name", "", d.gene_id or "", ""))
            continue
        if _cds_length(d) < cfg.min_cds_bp:
            decisions.append(Decision("denovo_skipped_cds_length", "", d.gene_id or "", f"cds<{cfg.min_cds_bp}"))
            continue
        candidates.append(d)

    liftoff_index = IntervalIndex((g.seqid, g.start, g.end, g) for g in liftoff_genes)
    removed: set = set()          # id(g) of Liftoff genes that were replaced
    placed: set = set()           # id(d) of de-novo genes placed via replace/add

    # ---- rule 2: replace (per-feature) -------------------------------------
    for d in candidates:
        did = id(d)
        overlappers = [
            l for l in liftoff_index.overlaps(d.seqid, d.start, d.end)
            if overlap_fraction(d.start, d.end, l.start, l.end) >= cfg.overlap_fraction
        ]
        replaced_any = False
        for l in overlappers:
            if id(l) in removed:
                continue
            if l.is_loc:
                removed.add(id(l))
                replaced_any = True
                d.evidence = "replace"
                decisions.append(Decision(
                    "replaced_loc", l.name, d.name,
                    f"{l.seqid}:{l.start}-{l.end} <- de-novo {d.gene_id}"))
            elif names_concordant(d.name, l.name):
                removed.add(id(l))
                replaced_any = True
                d.evidence = "replace"
                d.replace_concordant = True
                decisions.append(Decision(
                    "replaced_concordant", l.name, d.name,
                    f"{l.seqid}:{l.start}-{l.end} <- de-novo {d.gene_id}"))
            else:
                decisions.append(Decision(
                    "discordant_kept", l.name, d.name,
                    f"{l.seqid}:{l.start}-{l.end} kept (name differs)"))
        if replaced_any:
            placed.add(did)

    retained_liftoff = [g for g in liftoff_genes if id(g) not in removed]

    # ---- rule 3: add (novel, non-overlapping, non-duplicate) ----------------
    cds_index = IntervalIndex(
        (g.seqid, s, e, g) for g in retained_liftoff for (s, e) in g.cds_ranges)
    retained_base_names = {base_name(g.name) for g in retained_liftoff if g.name}

    for d in candidates:
        did = id(d)
        if did in placed:
            continue
        if not d.cds_ranges:
            decisions.append(Decision("denovo_skipped_no_cds", "", d.gene_id or "", ""))
            continue
        # strict CDS-overlap check against retained Liftoff genes
        overlaps_cds = False
        for s, e in d.cds_ranges:
            for l in cds_index.overlaps(d.seqid, s, e):
                if cds_overlap_fraction(d.cds_ranges, l.cds_ranges) > cfg.cds_overlap_fraction:
                    overlaps_cds = True
                    break
            if overlaps_cds:
                break
        if overlaps_cds:
            decisions.append(Decision("denovo_skipped_cds_overlap", "", d.gene_id or "", ""))
            continue
        if cfg.skip_duplicate_base_name and base_name(d.name) in retained_base_names:
            decisions.append(Decision("denovo_skipped_dup_base_name", "", d.gene_id or "", d.name))
            continue
        if cfg.bigwig is not None and _expression_mean(d, cfg) < cfg.min_expression:
            decisions.append(Decision("denovo_skipped_expression", "", d.gene_id or "", ""))
            continue
        placed.add(did)
        d.evidence = "de_novo"
        decisions.append(Decision(
            "added", "", d.name,
            f"{d.seqid}:{d.start}-{d.end} (de-novo {d.gene_id})"))

    # ---- assemble -----------------------------------------------------------
    final_genes: List[Gene] = list(retained_liftoff)
    for d in candidates:
        if id(d) in placed:
            d.set_name(d.name)          # normalise ID=gene-<name>;Name=<name>
            final_genes.append(d)
    assign_names(final_genes)
    for g in final_genes:
        g.normalize_children()  # rebuild Parent linkage (de-novo genes carry stale/missing Parent)

    return MergeResult(final_genes=final_genes, decisions=decisions)
