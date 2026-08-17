#!/usr/bin/env python3
"""
Merge HiCAT chr25 sub-chunk HTRM outputs into unified centromere annotation.

chr25 = 60,231,632 bp. Chunking: 5 parts (~12 Mb each) × 10 sub-chunks
(~1.2 Mb each, ~94k blocks) = 50 sub-chunks. 10 kb overlap.

Some sub-chunks were re-split into a/b halves for rescue (47k blocks each).
This script detects and handles both original and a/b rescue outputs.

Usage:
  python merge_chr25.py          # merge HTRM + decomposition TSVs
  python merge_chr25.py --tsv    # merge decomposition TSVs only

Output: HiCAT_genome/chr25_merged/
"""

import sys
from pathlib import Path
from collections import defaultdict

# ── Paths ───────────────────────────────────────────────────────────
BASE_DIR = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
                "/code/command-line-script/genome-annotation/HiCAT")
GENOME_DIR = BASE_DIR / "HiCAT_genome"
OUT_DIR = GENOME_DIR / "chr25_merged"

# ── Split parameters ────────────────────────────────────────────────
CHR25_TOTAL = 60_231_632
N_PARTS = 5
N_SUBS = 10
OVERLAP = 10_000
PART_LEN = CHR25_TOTAL // N_PARTS  # 12,046,326 (remainder in last part)


def build_coord_table():
    """Compute chr25 coordinate offsets for every sub-chunk region.

    Returns list of dicts. Each dict may have:
      - label, dir, offset, chr_start, chr_end, core_start, core_end
      - sub_total: sub-chunk FASTA length (for a/b split point)
      - is_rescue, half_a, half_b: for a/b rescue entries
    """
    rows = []
    for p in range(N_PARTS):
        part_num = p + 1  # 1-indexed for naming
        raw_part_start = p * PART_LEN
        raw_part_end = (p + 1) * PART_LEN if p < N_PARTS - 1 else CHR25_TOTAL
        part_start = max(0, raw_part_start - OVERLAP)
        part_end = min(CHR25_TOTAL, raw_part_end + OVERLAP)
        part_total = part_end - part_start

        SUB_LEN = part_total // N_SUBS

        for s in range(N_SUBS):
            sub_num = s + 1  # 1-indexed
            raw_sub_start = s * SUB_LEN
            raw_sub_end = (s + 1) * SUB_LEN if s < N_SUBS - 1 else part_total
            sub_start = max(0, raw_sub_start - OVERLAP)
            sub_end = min(part_total, raw_sub_end + OVERLAP)
            sub_total = sub_end - sub_start

            label = f"chr25_part{part_num:02d}_sub{sub_num:02d}"
            sub_dir = GENOME_DIR / label

            # chr25 genomic offset
            offset = part_start + sub_start

            # Core region (non-overlap) in chr25 coords
            # Clamp to part boundaries: sub-chunk cores should not extend into
            # part-level overlap regions (which are also covered by adjacent parts)
            core_chr_start = max(raw_part_start, part_start + raw_sub_start)
            core_chr_end = min(raw_part_end, part_start + raw_sub_end)

            entry = {
                "label": label,
                "dir": sub_dir,
                "offset": offset,
                "chr_start": part_start,
                "chr_end": part_end,
                "core_start": core_chr_start,
                "core_end": core_chr_end,
                "sub_total": sub_total,
            }

            # Determine if original completed or we need a/b rescue
            has_original = (sub_dir / "out" / "out_all_layer.xls").exists()
            dir_a = GENOME_DIR / f"{label}a"
            dir_b = GENOME_DIR / f"{label}b"
            has_a = (dir_a / "out" / "out_all_layer.xls").exists()
            has_b = (dir_b / "out" / "out_all_layer.xls").exists()

            if has_original:
                # Use original (preferred — ran on contiguous block)
                entry["source"] = "original"
                rows.append(entry)
            elif has_a and has_b:
                # Use a/b rescue pair
                half = sub_total // 2
                # a: covers sub-chunk coords [0, half+OVERLAP)
                # b: covers sub-chunk coords [half-OVERLAP, sub_total)
                entry_a = {
                    "label": f"{label}a",
                    "dir": dir_a,
                    "offset": offset,  # a starts at sub-chunk start
                    "chr_start": part_start,
                    "chr_end": part_end,
                    "core_start": core_chr_start,
                    "core_end": min(core_chr_end, offset + half - OVERLAP),
                    "source": "rescue_a",
                    "parent_label": label,
                }
                entry_b = {
                    "label": f"{label}b",
                    "dir": dir_b,
                    "offset": offset + half - OVERLAP,  # b starts at half-OVERLAP in sub-chunk
                    "chr_start": part_start,
                    "chr_end": part_end,
                    "core_start": max(core_chr_start, offset + half),
                    "core_end": core_chr_end,
                    "source": "rescue_b",
                    "parent_label": label,
                }
                rows.append(entry_a)
                rows.append(entry_b)
            else:
                print(f"  WARNING: {label}: no original and incomplete a/b — SKIPPED")

    return rows


# ── HTRM output readers ─────────────────────────────────────────────

def read_all_layer(path):
    """Read out_all_layer.xls → list of dicts (local coords).

    Columns: HOR_name  start  end  repeat_number  pattern  type
    """
    rows = []
    if not path.exists():
        return rows
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            rows.append({
                "hor": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "repeat": int(parts[3]),
                "pattern": parts[4],
                "type": parts[5] if len(parts) > 5 else "top",
            })
    return rows


def read_fasta(path):
    """Read FASTA → list of (header, sequence)."""
    entries = []
    if not path.exists():
        return entries
    with open(path) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            entries.append((header, "".join(seq_lines)))
    return entries


# ── Overlap resolution ──────────────────────────────────────────────

def resolve_overlaps(hors):
    """Deduplicate overlapping HOR calls in boundary regions.

    Sorts by chr25 start. When two HORs overlap:
      - Same HOR name → merge (extend end, keep max repeat)
      - Different name → keep the one with larger repeat_number
    """
    if not hors:
        return hors

    hors.sort(key=lambda h: (h["chr25_start"], h["chr25_end"]))
    resolved = []
    for hor in hors:
        if not resolved:
            resolved.append(hor)
            continue
        prev = resolved[-1]
        if hor["chr25_start"] <= prev["chr25_end"]:
            if hor["hor"] == prev["hor"]:
                prev["chr25_end"] = max(prev["chr25_end"], hor["chr25_end"])
                prev["repeat"] = max(prev["repeat"], hor["repeat"])
                prev["source"] = prev.get("source", "") + "+" + hor.get("source", "")
            else:
                if hor["repeat"] > prev["repeat"]:
                    hor["chr25_start"] = prev["chr25_start"]
                    resolved[-1] = hor
        else:
            resolved.append(hor)
    return resolved


# ── Decomposition TSV merging ───────────────────────────────────────

VARIANTS = [
    "final_decomposition.tsv",
    "final_decomposition_raw.tsv",
    "final_decomposition_alt.tsv",
]


def merge_decomposition_tsvs(coord_table):
    """Merge final_decomposition*.tsv from all sub-chunks.

    Filters monomers by core region (no overlap) to avoid duplicates.
    Offsets local coordinates → chr25 coordinates.
    """
    for variant in VARIANTS:
        out_path = OUT_DIR / variant
        all_rows = []
        total_kept = 0
        total_skipped = 0

        for cinfo in coord_table:
            tsv_path = cinfo["dir"] / variant
            if not tsv_path.exists():
                continue

            # Core region in chunk-local coords (0-based)
            core_local_start = cinfo["core_start"] - cinfo["offset"]
            core_local_end = cinfo["core_end"] - cinfo["offset"]

            kept = 0
            skipped = 0
            with open(tsv_path) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    cols = line.split("\t")
                    if len(cols) < 4:
                        continue
                    local_start = int(cols[2])
                    chr_start_0b = cinfo["offset"] + local_start
                    local_end = int(cols[3])
                    chr_end_0b = cinfo["offset"] + local_end

                    if core_local_start <= local_start < core_local_end:
                        cols[0] = "chr25"
                        cols[2] = str(chr_start_0b)
                        cols[3] = str(chr_end_0b)
                        all_rows.append("\t".join(cols))
                        kept += 1
                    else:
                        skipped += 1

            source_tag = cinfo.get("source", "?")
            print(f"    {cinfo['label']} [{source_tag}]: {variant} kept={kept:,} skipped={skipped:,}")
            total_kept += kept
            total_skipped += skipped

        if all_rows:
            with open(out_path, "w") as f:
                for row in all_rows:
                    f.write(row + "\n")
            print(f"  → {variant}: {total_kept:,} rows ({total_skipped:,} overlap removed)")
        else:
            print(f"  → {variant}: no data (skipped)")


# ── Main merge ──────────────────────────────────────────────────────

def merge_htrm(coord_table):
    """Merge HTRM outputs from all sub-chunks (original + a/b rescue)."""
    all_hors = []
    chunk_contrib = {}

    for cinfo in coord_table:
        layer_path = cinfo["dir"] / "out" / "out_all_layer.xls"
        if not layer_path.exists():
            print(f"  SKIP {cinfo['label']}: out/out_all_layer.xls not found")
            continue

        local_hors = read_all_layer(layer_path)
        for h in local_hors:
            h["chr25_start"] = cinfo["offset"] + h["start"]
            h["chr25_end"] = cinfo["offset"] + h["end"]
            h["source"] = cinfo["label"]
            del h["start"], h["end"]

        all_hors.extend(local_hors)
        chunk_contrib[cinfo["label"]] = len(local_hors)
        source_tag = cinfo.get("source", "?")
        print(f"  LOAD {cinfo['label']} [{source_tag}]: {len(local_hors)} HORs")

    if not all_hors:
        print("\nNo completed sub-chunks found. Run HiCAT jobs first.")
        return False

    print(f"\n  Total collected: {len(all_hors)} HORs")
    resolved = resolve_overlaps(all_hors)
    print(f"  After dedup:     {len(resolved)} HORs "
          f"({len(all_hors) - len(resolved)} duplicates removed)")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # out_all_layer.xls (0-based chr25 coords)
    with open(OUT_DIR / "out_all_layer.xls", "w") as f:
        for h in resolved:
            f.write(f"{h['hor']}\t{h['chr25_start']}\t{h['chr25_end']}\t"
                    f"{h['repeat']}\t{h['pattern']}\t{h['type']}\n")
    print(f"  Wrote: out_all_layer.xls ({len(resolved)} HORs)")

    # out_top_layer.xls (type == "top" only)
    top_count = 0
    with open(OUT_DIR / "out_top_layer.xls", "w") as f:
        for h in resolved:
            if h["type"] == "top":
                f.write(f"{h['hor']}\t{h['chr25_start']}\t{h['chr25_end']}\t"
                        f"{h['repeat']}\t{h['pattern']}\n")
                top_count += 1
    print(f"  Wrote: out_top_layer.xls ({top_count} top-layer HORs)")

    # HOR FASTA files
    for fa_name in ["out_hor.raw.fa", "out_hor.normal.fa"]:
        merged, seen = [], set()
        for cinfo in coord_table:
            fa_path = cinfo["dir"] / "out" / fa_name
            if not fa_path.exists():
                continue
            for header, seq in read_fasta(fa_path):
                new_hdr = f"{header} chunk={cinfo['label']}"
                if new_hdr not in seen:
                    merged.append((new_hdr, seq))
                    seen.add(new_hdr)
        if merged:
            with open(OUT_DIR / fa_name, "w") as f:
                for hdr, seq in merged:
                    f.write(f">{hdr}\n")
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i+80] + "\n")
            print(f"  Wrote: {fa_name} ({len(merged)} seqs)")

    # hor.repeatnumber.xls (sum across chunks)
    repeat_counts = defaultdict(int)
    for cinfo in coord_table:
        rn_path = cinfo["dir"] / "out" / "hor.repeatnumber.xls"
        if not rn_path.exists():
            continue
        with open(rn_path) as f:
            next(f)  # skip header
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    repeat_counts[parts[0]] += int(parts[1])
    if repeat_counts:
        with open(OUT_DIR / "hor.repeatnumber.xls", "w") as f:
            f.write("HORs\tRepeatNumber\n")
            for hor in sorted(repeat_counts.keys()):
                f.write(f"{hor}\t{repeat_counts[hor]}\n")
        print(f"  Wrote: hor.repeatnumber.xls ({len(repeat_counts)} HOR types)")

    # sort_statistics.xls (average across contributing chunks)
    stats_layers = defaultdict(list)
    for cinfo in coord_table:
        ss_path = cinfo["dir"] / "out" / "sort_statistics.xls"
        if not ss_path.exists():
            continue
        with open(ss_path) as f:
            next(f)  # skip header
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    stats_layers[int(parts[2])].append({
                        "cov": float(parts[0]),
                        "max_cov": float(parts[1]),
                    })
    if stats_layers:
        with open(OUT_DIR / "sort_statistics.xls", "w") as f:
            f.write("cov_rate\tmax_cov_rate\tsimilarity\n")
            for k in sorted(stats_layers.keys()):
                vals = stats_layers[k]
                f.write(f"{sum(v['cov'] for v in vals)/len(vals)}\t"
                        f"{sum(v['max_cov'] for v in vals)/len(vals)}\t{k}\n")
        print(f"  Wrote: sort_statistics.xls ({len(stats_layers)} layers)")

    # Provenance README
    with open(OUT_DIR / "README", "w") as f:
        f.write("Merged HiCAT chr25 centromere annotation\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Chr25 length: {CHR25_TOTAL:,} bp\n")
        f.write(f"Parts: {N_PARTS} × ~{PART_LEN:,} bp + {OVERLAP:,} bp overlap\n")
        f.write(f"Sub-chunks per part: {N_SUBS}\n")
        f.write(f"Total regions: {N_PARTS * N_SUBS}\n\n")
        f.write(f"Centromere location (centroAnno): 6,607,000-26,270,637\n\n")

        # Count by source type
        n_original = sum(1 for c in coord_table if c.get("source") == "original")
        n_rescue = sum(1 for c in coord_table if c.get("source") in ("rescue_a", "rescue_b"))
        f.write(f"Source: {n_original} original sub-chunks + "
                f"{n_rescue} rescue halves ({(n_rescue)//2} sub-chunks via rescue)\n\n")

        f.write("Contributions:\n")
        for label, count in sorted(chunk_contrib.items()):
            f.write(f"  {label}: {count} HORs\n")
        f.write(f"\nMerged: {len(resolved)} HORs "
                f"({len(all_hors) - len(resolved)} duplicates removed)\n")
        f.write(f"\nMerge date: {__import__('datetime').datetime.now()}\n")
    print(f"  Wrote: README")

    # Summary
    print(f"\n{'='*60}")
    print(f"Merge complete → {OUT_DIR}")
    print(f"{'='*60}")
    if repeat_counts:
        for hor, count in sorted(repeat_counts.items(), key=lambda x: -x[1])[:5]:
            print(f"  {hor}: {count:,} repeats")
    return True


# ── CLI ─────────────────────────────────────────────────────────────

def main():
    tsv_only = "--tsv" in sys.argv

    print("=" * 60)
    print("HiCAT chr25 merger")
    print("=" * 60)
    print(f"\nChr25: {CHR25_TOTAL:,} bp")
    print(f"Parts: {N_PARTS} × ~{PART_LEN:,} bp")
    print(f"Sub-chunks: {N_SUBS} per part = {N_PARTS * N_SUBS} total")
    print(f"Overlap: {OVERLAP:,} bp\n")

    coord_table = build_coord_table()
    print(f"Coord table: {len(coord_table)} entries\n")

    # Show status
    n_original = sum(1 for c in coord_table if c.get("source") == "original")
    n_rescue_a = sum(1 for c in coord_table if c.get("source") == "rescue_a")
    n_rescue_b = sum(1 for c in coord_table if c.get("source") == "rescue_b")
    print(f"  Original sub-chunks: {n_original}")
    print(f"  Rescue halves (a+b): {n_rescue_a + n_rescue_b} (from {(n_rescue_a + n_rescue_b)//2} sub-chunks)")
    print()

    for c in coord_table:
        has_out = (c["dir"] / "out" / "out_all_layer.xls").exists()
        marker = " [READY]" if has_out else " [missing]"
        tag = c.get("source", "?")
        print(f"  {c['label']:35s} [{tag:10s}]  chr25:{c['offset']:>10,}  "
              f"core:{c['core_start']:>10,}-{c['core_end']:>10,}{marker}")

    if tsv_only:
        print(f"\n{'='*60}")
        print("Merging decomposition TSVs only")
        print(f"{'='*60}")
        merge_decomposition_tsvs(coord_table)
    else:
        print(f"\n{'='*60}")
        print("Merging HTRM outputs")
        print(f"{'='*60}")
        ok = merge_htrm(coord_table)
        if ok:
            print(f"\n{'='*60}")
            print("Also merging decomposition TSVs")
            print(f"{'='*60}")
            merge_decomposition_tsvs(coord_table)

    print(f"\nDone. Output: {OUT_DIR}/")


if __name__ == "__main__":
    main()
