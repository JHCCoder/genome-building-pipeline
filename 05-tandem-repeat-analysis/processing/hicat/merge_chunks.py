#!/usr/bin/env python3
"""
Merge HiCAT sub-chunk HTRM outputs into a unified chr4 centromere annotation.

Each sub-chunk runs HiCAT independently with local coordinates. This script:
  1. Computes the chr4 genomic offset for every sub-chunk
  2. Offsets all HOR coordinates to the chr4 coordinate system
  3. Resolves overlapping HOR calls in the 10 kb boundary regions
  4. Concatenates HOR FASTA sequences with unique identifiers
  5. Sums repeat-number statistics across chunks
  6. Regenerates sort_statistics.xls from the merged all_layer output
  7. Also merges decomposition TSVs (final_decomposition*.tsv)

Usage:
  python merge_chunks.py          # merge chr4 HTRM output
  python merge_chunks.py --tsv    # merge decomposition TSVs only

Output written to: HiCAT_genome/chr4_merged/

NOTE: Run this AFTER all sub-chunk jobs have completed (out/ directory present).
      Incomplete chunks are silently skipped.
"""

import sys
from pathlib import Path
from collections import defaultdict

# ── Paths ───────────────────────────────────────────────────────────
BASE_DIR = Path("/tscc/projects/ps-renlab2/jhc103/degu-genome-assembly-proj"
                "/code/command-line-script/genome-annotation/HiCAT")
GENOME_DIR = BASE_DIR / "HiCAT_genome"
OUT_DIR = GENOME_DIR / "chr4_merged"

# ── Split parameters (must match split_chr4.py & resplit_failed_chunks.py) ──
CHR4_TOTAL = 151_562_733
N_CHUNKS = 10
OVERLAP = 10_000
CHUNK_LEN = CHR4_TOTAL // N_CHUNKS  # 15,156,273
N_SUB = 2  # chr4 chunks split into 2 sub-chunks (except part09)


def build_coord_table():
    """Compute chr4 coordinate offsets for every sub-chunk.

    Returns list of dicts with:
      label, dir, offset (local → chr4 pos), chunk_start, chunk_end,
      core_start, core_end
    """
    rows = []
    for i in range(N_CHUNKS):
        part_num = i + 1
        raw_start = i * CHUNK_LEN
        raw_end = (i + 1) * CHUNK_LEN if i < N_CHUNKS - 1 else CHR4_TOTAL
        chunk_start = max(0, raw_start - OVERLAP)
        chunk_end = min(CHR4_TOTAL, raw_end + OVERLAP)
        chunk_total = chunk_end - chunk_start

        # chr4_part09 (centromere) was split into 8 sub-chunks (rescue round)
        # Sub-chunk fasta headers give chunk-local coords, e.g.:
        #   >chr4_part09_sub01 chr4_part09 chunk:1-1907034
        #   >chr4_part09_sub02 chr4_part09 chunk:1887035-3804068
        # offset = chunk_start + sub_local_start - 1
        if part_num == 9:
            part09_subs = [
                ("01", 1,          1907034),
                ("02", 1887035,    3804068),
                ("03", 3784069,    5701102),
                ("04", 5681103,    7598136),
                ("05", 7578137,    9495170),
                ("06", 9475171,   11392204),
                ("07", 11372205,  13289238),
                ("08", 13269239,  15176273),
            ]
            for sub_id, local_start, local_end in part09_subs:
                sub_offset = chunk_start + local_start - 1
                rows.append({
                    "label": f"chr4_part{part_num:02d}_sub{sub_id}",
                    "dir": GENOME_DIR / f"chr4_part{part_num:02d}_sub{sub_id}",
                    "chr_start": chunk_start,
                    "chr_end": chunk_end,
                    "core_start": raw_start,
                    "core_end": raw_end,
                    "offset": sub_offset,
                    "is_full_chunk": False,
                })
        else:
            sub_len = chunk_total // N_SUB
            for j in range(N_SUB):
                raw_sub_start = j * sub_len
                raw_sub_end = (j + 1) * sub_len if j < N_SUB - 1 else chunk_total
                sub_start = max(0, raw_sub_start - OVERLAP)
                sub_end = min(chunk_total, raw_sub_end + OVERLAP)

                rows.append({
                    "label": f"chr4_part{part_num:02d}_sub{j+1:02d}",
                    "dir": GENOME_DIR / f"chr4_part{part_num:02d}_sub{j+1:02d}",
                    "chr_start": chunk_start,
                    "chr_end": chunk_end,
                    "core_start": raw_start,
                    "core_end": raw_end,
                    "offset": chunk_start + sub_start,
                    "is_full_chunk": False,
                })
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

    Sorts by chr4 start. When two HORs overlap:
      - Same HOR name → merge (extend end, keep max repeat)
      - Different name → keep the one with larger repeat_number
    """
    if not hors:
        return hors

    hors.sort(key=lambda h: (h["chr4_start"], h["chr4_end"]))
    resolved = []
    for hor in hors:
        if not resolved:
            resolved.append(hor)
            continue
        prev = resolved[-1]
        if hor["chr4_start"] <= prev["chr4_end"]:
            if hor["hor"] == prev["hor"]:
                prev["chr4_end"] = max(prev["chr4_end"], hor["chr4_end"])
                prev["repeat"] = max(prev["repeat"], hor["repeat"])
                prev["source"] = prev.get("source", "") + "+" + hor.get("source", "")
            else:
                if hor["repeat"] > prev["repeat"]:
                    hor["chr4_start"] = prev["chr4_start"]
                    resolved[-1] = hor
        else:
            resolved.append(hor)
    return resolved


# ── Decomposition TSV merging (with overlap filtering) ──────────────

VARIANTS = [
    "final_decomposition.tsv",
    "final_decomposition_raw.tsv",
    "final_decomposition_alt.tsv",
]


def merge_decomposition_tsvs(coord_table):
    """Merge final_decomposition*.tsv from all sub-chunks.

    Filters out monomers in overlap regions (keeps only core-region monomers).
    Offsets local coordinates → chr4 coordinates.
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

            # core region in chunk-local coords
            core_local_start = cinfo["core_start"] - cinfo["chr_start"]
            core_local_end = cinfo["core_end"] - cinfo["chr_start"]

            kept = 0
            skipped = 0
            with open(tsv_path) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    cols = line.split("\t")
                    local_start = int(cols[2])
                    chr_start_0b = cinfo["offset"] + local_start
                    chr_start_1b = chr_start_0b + 1

                    if core_local_start <= local_start <= core_local_end:
                        cols[0] = "chr4"
                        cols[2] = str(chr_start_0b)
                        local_end = int(cols[3])
                        chr_end_0b = cinfo["offset"] + local_end
                        cols[3] = str(chr_end_0b)
                        all_rows.append("\t".join(cols))
                        kept += 1
                    else:
                        skipped += 1

            print(f"    {cinfo['label']}: {variant} kept={kept:,} skipped={skipped:,}")
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
    """Merge HTRM outputs from all sub-chunks."""
    all_hors = []
    chunk_contrib = {}

    for cinfo in coord_table:
        layer_path = cinfo["dir"] / "out" / "out_all_layer.xls"
        if not layer_path.exists():
            print(f"  SKIP {cinfo['label']}: out/out_all_layer.xls not found")
            continue

        local_hors = read_all_layer(layer_path)
        for h in local_hors:
            h["chr4_start"] = cinfo["offset"] + h["start"]
            h["chr4_end"] = cinfo["offset"] + h["end"]
            h["source"] = cinfo["label"]
            del h["start"], h["end"]

        all_hors.extend(local_hors)
        chunk_contrib[cinfo["label"]] = len(local_hors)
        print(f"  LOAD {cinfo['label']}: {len(local_hors)} HORs")

    if not all_hors:
        print("\nNo completed sub-chunks found. Run HiCAT jobs first.")
        return False

    print(f"\n  Total collected: {len(all_hors)} HORs")
    resolved = resolve_overlaps(all_hors)
    print(f"  After dedup:     {len(resolved)} HORs "
          f"({len(all_hors) - len(resolved)} duplicates removed)")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # out_all_layer.xls
    with open(OUT_DIR / "out_all_layer.xls", "w") as f:
        for h in resolved:
            f.write(f"{h['hor']}\t{h['chr4_start']}\t{h['chr4_end']}\t"
                    f"{h['repeat']}\t{h['pattern']}\t{h['type']}\n")
    print(f"  Wrote: out_all_layer.xls ({len(resolved)} HORs)")

    # out_top_layer.xls (type == "top" only)
    top_count = 0
    with open(OUT_DIR / "out_top_layer.xls", "w") as f:
        for h in resolved:
            if h["type"] == "top":
                f.write(f"{h['hor']}\t{h['chr4_start']}\t{h['chr4_end']}\t"
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

    # Provenance
    with open(OUT_DIR / "README", "w") as f:
        f.write("Merged HiCAT chr4 centromere annotation\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Chr4 length: {CHR4_TOTAL:,} bp\n")
        f.write(f"Chunks: {N_CHUNKS} × ~{CHUNK_LEN:,} bp + {OVERLAP:,} bp overlap\n")
        f.write(f"Sub-chunks per chunk: {N_SUB} (except part09, full chunk)\n\n")
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
    print("HiCAT chr4 merger")
    print("=" * 60)

    coord_table = build_coord_table()
    print(f"\nCoord table: {len(coord_table)} entries\n")
    for c in coord_table:
        has_out = (c["dir"] / "out" / "out_all_layer.xls").exists()
        marker = " [READY]" if has_out else " [pending]"
        print(f"  {c['label']:30s}  chr4:{c['chr_start']:,}-{c['chr_end']:,}"
              f"  offset={c['offset']:>9,}{marker}")

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
