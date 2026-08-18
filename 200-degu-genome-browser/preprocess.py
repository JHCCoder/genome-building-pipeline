#!/usr/bin/env python3
"""
Pre-compute static data files for the Apache-only genome browser.
Run this ONCE on silencer (where FASTA, BigWigs, and Hi-C live).
Output goes into static/ — copy that folder alongside index.html.

Requires: conda activate browser  (has pyfaidx, pyBigWig, h5py, hdf5plugin, numpy)
Usage:    python3 preprocess.py
"""

import os, sys, json, struct, math, bisect
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
STATIC = os.path.join(HERE, "static")

# ---- Paths — adjust DEPLOY to where your data files live ----
DEPLOY = os.path.join(HERE, "..", "degu-genome-browser")
FASTA_PATH = os.path.join(DEPLOY, "assembly_final.sorted.headerRenamed.chrAssigned.mito.fasta")
BIGWIG_PATHS = {
    "cenpa": os.path.join(DEPLOY, "XG_150.all.bw"),
    "h3k27ac": os.path.join(DEPLOY, "XG_152.all.bw"),
}
HIC_PATH = os.path.join(DEPLOY, "403_ear_deep_file_test_hifi_041425_10kb.h5")
HOR_PATH = os.path.join(DEPLOY, "HORs.bed")
HICAT_PATH = os.path.join(DEPLOY, "hicat.bed")  # HiCAT chr4:124-134Mb centromere HOR units
SEGDUP_PATH = os.path.join(DEPLOY, "segdup_output_mod.bedpe")
TANDEM_DIR = os.path.join(DEPLOY, "tandem_data")

# Find actual FASTA (resolve symlink)
for p in [FASTA_PATH]:
    if not os.path.exists(p):
        # try resolve symlink
        real = os.path.realpath(p)
        if os.path.exists(real):
            FASTA_PATH = real

print("=" * 60)
print("Static data pre-processor for genome browser")
print("=" * 60)

# ============================================================================
# 1. Sequence — plain text per chromosome (1 byte per base)
# ============================================================================
print("\n[1/6] Sequence files...")
import pyfaidx

genome = pyfaidx.Fasta(FASTA_PATH)
seq_dir = os.path.join(STATIC, "sequence")
os.makedirs(seq_dir, exist_ok=True)

for chrom in sorted(genome.keys()):
    out_path = os.path.join(seq_dir, f"{chrom}.seq")
    if os.path.exists(out_path):
        print(f"  {chrom}: already exists, skipping")
        continue
    seq = str(genome[chrom][:]).upper()
    with open(out_path, "w") as f:
        f.write(seq)
    size_mb = len(seq) / 1e6
    print(f"  {chrom}: {len(seq):,} bp ({size_mb:.1f} MB)")

print(f"  Done: {len(genome.keys())} chromosomes → {seq_dir}/")

# ============================================================================
# 2. BigWig — pre-computed binned means at 100bp resolution (binary float32)
# ============================================================================
print("\n[2/6] BigWig signal tracks...")
import pyBigWig
import numpy as np

BIN_SIZE = 100  # bp per bin — client can downsample from here

for track_name, bw_path in BIGWIG_PATHS.items():
    if not os.path.exists(bw_path):
        print(f"  {track_name}: file not found at {bw_path}, skipping")
        continue

    bw = pyBigWig.open(bw_path)
    out_dir = os.path.join(STATIC, "bigwig", track_name)
    os.makedirs(out_dir, exist_ok=True)

    for chrom in sorted(genome.keys()):
        out_path = os.path.join(out_dir, f"{chrom}.bin")
        if os.path.exists(out_path):
            print(f"  {track_name}/{chrom}: already exists, skipping")
            continue

        chrom_len = len(genome[chrom])
        n_bins = (chrom_len + BIN_SIZE - 1) // BIN_SIZE

        try:
            stats = bw.stats(chrom, 0, chrom_len, nBins=n_bins,
                             type="mean", exact=True)
        except Exception as e:
            print(f"  {track_name}/{chrom}: error — {e}, filling with zeros")
            stats = [0.0] * n_bins

        # Replace None/NaN with -1 sentinel (no data)
        values = np.array([v if (v is not None and v == v) else -1.0
                          for v in stats], dtype=np.float32)

        with open(out_path, "wb") as f:
            f.write(struct.pack("<I", BIN_SIZE))   # header: bin size
            f.write(struct.pack("<I", chrom_len))   # header: chrom length
            f.write(values.tobytes())

        print(f"  {track_name}/{chrom}: {n_bins:,} bins, "
              f"{values.nbytes / 1024:.0f} KB")

    bw.close()

print(f"  Done: BigWig tracks → {STATIC}/bigwig/")

# ============================================================================
# 3. Hi-C — pre-computed sparse contacts at 10kb resolution (per-chrom JSON)
# ============================================================================
print("\n[3/6] Hi-C contact matrices...")
import hdf5plugin
import h5py

if not os.path.exists(HIC_PATH):
    print(f"  Hi-C file not found at {HIC_PATH}, skipping")
else:
    hic_dir = os.path.join(STATIC, "hic")
    os.makedirs(hic_dir, exist_ok=True)

    h5 = h5py.File(HIC_PATH, "r")
    data = h5["matrix/data"][:]
    indices = h5["matrix/indices"][:]
    indptr = h5["matrix/indptr"][:]
    chr_raw = h5["intervals/chr_list"][:]
    chr_names = [c.decode("utf-8") for c in chr_raw]
    starts = h5["intervals/start_list"][:]
    ends = h5["intervals/end_list"][:]

    # Build {chrom: (first_idx, last_idx+1)}
    chrom_ranges = {}
    n = len(chr_names)
    i = 0
    while i < n:
        chrom = chr_names[i]
        first = i
        while i < n and chr_names[i] == chrom:
            i += 1
        chrom_ranges[chrom] = (first, i)

    for chrom in sorted(chrom_ranges.keys()):
        out_path = os.path.join(hic_dir, f"{chrom}.json")
        if os.path.exists(out_path):
            print(f"  hic/{chrom}: already exists, skipping")
            continue

        first, last = chrom_ranges[chrom]
        n_bins = last - first
        bin_size = int(ends[first] - starts[first])

        # Extract sparse upper-triangle contacts
        contacts = []
        for row in range(first, last):
            local_row = row - first
            row_start = int(indptr[row])
            row_end = int(indptr[row + 1])
            for p in range(row_start, row_end):
                col = int(indices[p])
                if col < first:
                    continue
                if col >= last:
                    break
                if col >= row:  # upper triangle
                    val = int(data[p])
                    if val > 0:
                        contacts.append([local_row, col - first, val])

        bin_starts = [int(starts[i]) for i in range(first, last)]
        bin_ends = [int(ends[i]) for i in range(first, last)]

        out = {
            "chrom": chrom,
            "bin_size": bin_size,
            "n_bins": n_bins,
            "bins": [[bin_starts[i], bin_ends[i]] for i in range(n_bins)],
            "contacts": contacts,
        }

        with open(out_path, "w") as f:
            json.dump(out, f)

        size_kb = os.path.getsize(out_path) / 1024
        print(f"  hic/{chrom}: {n_bins} bins, {len(contacts):,} contacts, "
              f"{size_kb:.0f} KB")

    h5.close()
    print(f"  Done: Hi-C → {hic_dir}/")

# ============================================================================
# 4. HOR — convert BED to per-chromosome JSON
# ============================================================================
print("\n[4/6] Higher-order repeats...")
hor_dir = os.path.join(STATIC, "hor")
os.makedirs(hor_dir, exist_ok=True)

if not os.path.exists(HOR_PATH):
    print(f"  HOR file not found at {HOR_PATH}, skipping")
else:
    raw = defaultdict(list)
    with open(HOR_PATH) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            chrom = parts[0]
            start_0b = int(parts[2])
            end_0b = int(parts[3])
            monomers = int(parts[4])
            span = int(parts[6])
            name = parts[1]
            raw[chrom].append([start_0b + 1, end_0b, name, monomers, span])

    total = 0
    for chrom in sorted(raw):
        raw[chrom].sort(key=lambda x: x[0])
        out_path = os.path.join(hor_dir, f"{chrom}.json")
        with open(out_path, "w") as f:
            json.dump(raw[chrom], f)
        total += len(raw[chrom])
        print(f"  hor/{chrom}: {len(raw[chrom])} entries")

    print(f"  Done: {len(raw)} chromosomes, {total} entries → {hor_dir}/")

# ----------------------------------------------------------------------------
# 4b. HiCAT — centromere HOR annotation (chr4:124-134 Mb only)
# BED: chrom, start(0b), end(0b), RnLn, monomer_pattern, strand
# ----------------------------------------------------------------------------
print("\n[4b] HiCAT centromere HOR annotation...")
hicat_dir = os.path.join(STATIC, "hicat")
os.makedirs(hicat_dir, exist_ok=True)

if not os.path.exists(HICAT_PATH):
    print(f"  HiCAT file not found at {HICAT_PATH}, skipping")
else:
    raw = defaultdict(list)
    with open(HICAT_PATH) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start_0b = int(parts[1])
            end_0b = int(parts[2])
            rlnln = parts[3]
            pattern = parts[4]
            strand = parts[5]
            raw[chrom].append([start_0b + 1, end_0b, rlnln, pattern, strand])

    total = 0
    for chrom in sorted(raw):
        raw[chrom].sort(key=lambda x: x[0])
        out_path = os.path.join(hicat_dir, f"{chrom}.json")
        with open(out_path, "w") as f:
            json.dump(raw[chrom], f)
        total += len(raw[chrom])
        print(f"  hicat/{chrom}: {len(raw[chrom])} entries")

    print(f"  Done: {len(raw)} chromosomes, {total} entries → {hicat_dir}/")

# ============================================================================
# 5. Segmental duplications — convert BEDPE to per-chromosome JSON
# ============================================================================
print("\n[5/6] Segmental duplications...")
segdup_dir = os.path.join(STATIC, "segdup")
os.makedirs(segdup_dir, exist_ok=True)

if not os.path.exists(SEGDUP_PATH):
    print(f"  Segdup file not found at {SEGDUP_PATH}, skipping")
else:
    raw = defaultdict(list)
    with open(SEGDUP_PATH) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            c1, s1, e1 = parts[0], int(parts[1]), int(parts[2])
            c2, s2, e2 = parts[3], int(parts[4]), int(parts[5])
            s1 += 1  # 0-based → 1-based closed
            s2 += 1
            raw[c1].append([s1, e1, c2, s2, e2])
            raw[c2].append([s2, e2, c1, s1, e1])

    total = 0
    for chrom in sorted(raw):
        raw[chrom].sort(key=lambda x: x[0])
        out_path = os.path.join(segdup_dir, f"{chrom}.json")
        with open(out_path, "w") as f:
            json.dump(raw[chrom], f)
        total += len(raw[chrom])
        print(f"  segdup/{chrom}: {len(raw[chrom])} entries")

    print(f"  Done: {len(raw)} chromosomes, {total} entries → {segdup_dir}/")

# ============================================================================
# 6. Copy static directories from deploy
# ============================================================================
print("\n[6/6] Copying features & tandem data...")
import shutil

# features_genes/
src = os.path.join(DEPLOY, "features_genes")
dst = os.path.join(STATIC, "features_genes")
if os.path.exists(src) and not os.path.exists(dst):
    shutil.copytree(src, dst)
    print(f"  features_genes/ copied ({len(os.listdir(dst))} files)")

# features/
src = os.path.join(DEPLOY, "features")
dst = os.path.join(STATIC, "features")
if os.path.exists(src) and not os.path.exists(dst):
    shutil.copytree(src, dst)
    print(f"  features/ copied ({len(os.listdir(dst))} files)")

# tandem_data/
src = os.path.join(DEPLOY, "tandem_data")
dst = os.path.join(STATIC, "tandem")
if os.path.exists(src) and not os.path.exists(dst):
    shutil.copytree(src, dst)
    print(f"  tandem/ copied ({len(os.listdir(dst))} files)")

# structural_errors.json
src = os.path.join(DEPLOY, "structural_errors.json")
dst = os.path.join(STATIC, "structural_errors.json")
if os.path.exists(src) and not os.path.exists(dst):
    shutil.copy2(src, dst)
    print(f"  structural_errors.json copied")

# segdup_data.js (inline script)
src = os.path.join(DEPLOY, "segdup_data.js")
dst = os.path.join(HERE, "segdup_data.js")
if os.path.exists(src) and not os.path.exists(dst):
    shutil.copy2(src, dst)
    print(f"  segdup_data.js copied")

print("\n" + "=" * 60)
print("Preprocessing complete!")
print(f"Static data is in: {STATIC}/")
print(f"Copy the entire '{HERE}' folder to Apache.")
print("=" * 60)
