"""Gene-name utilities for the annotation merge.

Handles the conventions used by Liftoff output and Braker/functional-annotation
output:

* unannotated ("LOC") genes: names matching ``LOC\\d{9}`` (case-insensitive).
* predicted paralog copies: names carrying a trailing ``-(D?R?L\\d+)`` suffix
  (e.g. ``CLEC4A-rl1``, ``MAGEA10-DL1``).  Concordance is assessed on the
  *suffix-stripped* base name so that a de-novo copy of a gene is recognised as
  the same locus as the transferred gene it duplicates.
"""
from __future__ import annotations

import re

_LOC_RE = re.compile(r"^LOC\d{9}$", re.IGNORECASE)
# Trailing paralog-copy suffix: -L1 / -DL1 / -RL1 (any number of digits).
_SUFFIX_RE = re.compile(r"-(D?R?L\d+)$", re.IGNORECASE)
# The letters of a paralog suffix only (L / DL / RL), for type comparison.
_SUFFIX_TYPE_RE = re.compile(r"-(D?R?L)\d+$", re.IGNORECASE)
# BLAST/UniProt name attributes sometimes carry several names at once,
# e.g. "PCDHAC2 - PCDHA10 - PCDHA4" or "TRAPPC2B,TRAPPC2".
_MULTI_NAME_RE = re.compile(r"\s*,\s*|\s*;\s*|\s+-\s+")


def normalize(name: str | None) -> str:
    """Return an uppercase, stripped form of *name* for comparisons."""
    if not name:
        return ""
    return str(name).strip().upper()


def base_name(name: str | None) -> str:
    """Return the suffix-stripped, uppercased gene name.

    ``CLEC4A-rl1`` and ``CLEC4A`` both reduce to ``CLEC4A``.
    """
    n = normalize(name)
    return _SUFFIX_RE.sub("", n)


def is_loc(name: str | None) -> bool:
    """True if *name* is an unannotated ``LOC##########`` gene."""
    return bool(_LOC_RE.match(normalize(name)))


def has_suffix(name: str | None) -> bool:
    """True if *name* carries a paralog-copy suffix (``-L*``/``-DL*``/``-RL*``)."""
    return bool(_SUFFIX_RE.search(name or ""))


def suffix_type(name: str | None) -> str:
    """Return the paralog suffix *type* of *name* ('L', 'DL', 'RL') or ''."""
    m = _SUFFIX_TYPE_RE.search(name or "")
    return m.group(1).upper() if m else ""


def name_aliases(name: str | None) -> list:
    """Split a multi-name attribute (e.g. 'PCDHAC2 - PCDHA10 - PCDHA4') into
    its constituent gene names. A single name is returned unchanged."""
    return [x.strip() for x in _MULTI_NAME_RE.split(name or "") if x.strip()]


def base_name_aliases(name: str | None) -> set:
    """Set of suffix-stripped base names carried by a (possibly multi-)name."""
    return {base_name(x) for x in name_aliases(name)}


def names_concordant(a: str | None, b: str | None) -> bool:
    """True if *a* and *b* are the same gene (suffix-insensitive, case-insensitive).

    Multi-name attributes count as concordant if *any* of their names match,
    so a de-novo gene annotated 'PCDHAC2 - PCDHA10 - PCDHA4' is recognised as a
    PCDHA4 prediction."""
    return bool(base_name_aliases(a) & base_name_aliases(b))
