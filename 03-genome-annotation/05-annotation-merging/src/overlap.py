"""Interval-index helpers for finding gene / CDS overlaps efficiently."""
from __future__ import annotations

import bisect
from typing import Dict, Iterable, List, Tuple


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


class IntervalIndex:
    """Per-sequence sorted index of (start, end) intervals.

    Supports ``overlappers(query_start, query_end)`` returning every stored
    interval that overlaps the query, plus the stored payload.
    """

    def __init__(self, intervals: Iterable[Tuple[str, int, int, object]]):
        # group by seqid; keep (start, end, payload)
        self._by_seq: Dict[str, List[Tuple[int, int, object]]] = {}
        for seqid, start, end, payload in intervals:
            self._by_seq.setdefault(seqid, []).append((start, end, payload))
        self._starts: Dict[str, List[int]] = {}
        for seqid, items in self._by_seq.items():
            items.sort(key=lambda x: (x[0], x[1]))
            self._starts[seqid] = [x[0] for x in items]

    def overlaps(self, seqid: str, start: int, end: int):
        """Yield payloads of stored intervals overlapping [start, end]."""
        starts = self._starts.get(seqid)
        if not starts:
            return
        items = self._by_seq[seqid]
        lo = bisect.bisect_right(starts, end)  # first index with start > end
        for i in range(lo - 1, -1, -1):
            s, e, payload = items[i]
            if e < start:
                break
            if s <= end:  # overlap
                yield payload


def overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """Fraction of the *shorter* interval covered by the overlap."""
    ov = _overlap(a_start, a_end, b_start, b_end)
    if ov == 0:
        return 0.0
    denom = min(a_end - a_start + 1, b_end - b_start + 1)
    return ov / denom if denom else 0.0


def cds_overlap_fraction(cds_a: List[Tuple[int, int]],
                         cds_b: List[Tuple[int, int]]) -> float:
    """Fraction of the shorter CDS set's total length overlapped by the other."""
    total_a = sum(e - s + 1 for s, e in cds_a)
    total_b = sum(e - s + 1 for s, e in cds_b)
    if not total_a or not total_b:
        return 0.0
    overlap = 0
    for s1, e1 in cds_a:
        for s2, e2 in cds_b:
            overlap += _overlap(s1, e1, s2, e2)
    denom = min(total_a, total_b)
    return overlap / denom if denom else 0.0
