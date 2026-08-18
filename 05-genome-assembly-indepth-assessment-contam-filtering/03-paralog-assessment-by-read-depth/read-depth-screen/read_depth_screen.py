#!/usr/bin/env python3
"""
Read-depth haplotig screening for parent–paralog gene pairs.

Estimates genome-wide modal coverage from non-repetitive 1-kb windows,
then calculates per-base coverage (unique MAPQ≥20 and permissive) across
gene body ± flanking regions for each parent and paralog.  Pairs are
classified as:

  likely_genuine_duplication — both flanks near diploid coverage C
  suspicious_haplotig          — both loci at C/2, flanks included
  recent_paralog_mapping_ambig — gene body low unique but permissive OK
  inconclusive                 — insufficient mappable flank, or conflicting

The key design rule: a locus is NOT called a haplotig solely because the
gene body has C/2 unique coverage — half-depth must also appear in
informative flanking sequence.

Usage:
    python3 read_depth_screen.py \\
        --bam aligned.bam \\
        --fasta assembly.fasta \\
        --paralog-families paralog_families.tsv \\
        --genome-lengths genome_length.txt \\
        --out-dir ./output \\
        --repeat-bed repeats.gff \\
        --threads 24 \\
        --flank-bp 25000
"""

import argparse
import os
import sys
import time
import subprocess
import tempfile
import shutil
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

# ============================================================
# Constants
# ============================================================
AUTO_CHROMS = [f'chr{i}' for i in range(1, 29)]
WINDOW_SIZE = 1000  # 1-kb windows for genome-wide coverage estimation

# Classification labels
CLASS_GENUINE = 'likely_genuine_duplication'
CLASS_HAPLOTIG = 'suspicious_haplotig'
CLASS_AMBIG = 'recent_paralog_mapping_ambiguity'
CLASS_INCONCLUSIVE = 'inconclusive'

# ============================================================
# Helpers
# ============================================================
def find_coverage_mode(depths, min_count=100):
    """
    Find the main coverage peak (mode) from an array of window-mean depths.
    Uses KDE to locate the highest peak.
    Falls back to histogram mode if KDE fails.
    """
    depths = depths[np.isfinite(depths)]
    depths = depths[depths > 0]
    if len(depths) < min_count:
        return np.nan

    # Clip extremes for robust mode detection (0.5–99.5 percentile)
    lo, hi = np.percentile(depths, [0.5, 99.5])
    clipped = depths[(depths >= lo) & (depths <= hi)]

    try:
        # Use KDE to find the main peak
        bw = max(0.5, (hi - lo) / 100)  # bandwidth
        kde = gaussian_kde(clipped, bw_method=bw / np.std(clipped))
        x_grid = np.linspace(lo, hi, 500)
        y_grid = kde(x_grid)

        # Find peaks
        peaks, props = find_peaks(y_grid, height=np.max(y_grid) * 0.1)
        if len(peaks) > 0:
            # Return the highest peak
            best_idx = peaks[np.argmax(props['peak_heights'])]
            return x_grid[best_idx]
    except Exception:
        pass

    # Fallback: histogram mode
    n_bins = max(20, int(np.sqrt(len(clipped))))
    counts, edges = np.histogram(clipped, bins=n_bins)
    mode_bin = np.argmax(counts)
    return (edges[mode_bin] + edges[mode_bin + 1]) / 2


def load_repeat_mask(repeat_path, min_contig_len=0):
    """
    Load a repeat annotation file (GFF or BED) and return a dict
    mapping chrom -> list of (start, end) repeat intervals.
    Handles both GFF (col 4) and BED (col 1) naming conventions.
    """
    repeats = defaultdict(list)
    if not repeat_path or not os.path.exists(repeat_path):
        return repeats

    with open(repeat_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            chrom = parts[0]
            try:
                start = int(parts[3]) if len(parts) >= 5 else int(parts[1])
                end = int(parts[4]) if len(parts) >= 5 else int(parts[2])
            except ValueError:
                continue
            repeats[chrom].append((start, end))

    # Sort and optionally merge intervals
    for chrom in repeats:
        repeats[chrom].sort()
    return repeats


def window_overlaps_repeat(chrom, w_start, w_end, repeats):
    """Check if a window overlaps any repeat interval."""
    if chrom not in repeats:
        return False
    intervals = repeats[chrom]
    # Binary search
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        if intervals[mid][1] <= w_start:
            lo = mid + 1
        else:
            hi = mid
    # Check the intervals starting at `lo`
    for i in range(lo, len(intervals)):
        r_start, r_end = intervals[i]
        if r_start >= w_end:
            break
        if r_end > w_start and r_start < w_end:
            return True
    return False


def build_windows(chrom_lengths, repeats=None, window_size=WINDOW_SIZE):
    """
    Build 1-kb windows across autosomes, optionally excluding repeat-overlapping windows.
    Returns DataFrame with columns: chrom, start, end, is_repeat_free
    """
    windows = []
    for chrom in AUTO_CHROMS:
        if chrom not in chrom_lengths:
            continue
        chr_len = chrom_lengths[chrom]
        for w_start in range(0, chr_len, window_size):
            w_end = min(w_start + window_size, chr_len)
            is_rpt = False
            if repeats:
                is_rpt = window_overlaps_repeat(chrom, w_start, w_end, repeats)
            windows.append({
                'chrom': chrom,
                'start': w_start,
                'end': w_end,
                'length': w_end - w_start,
                'overlaps_repeat': is_rpt,
            })
    df = pd.DataFrame(windows)
    print(f"  Built {len(df):,} windows ({len(df[~df['overlaps_repeat']]):,} non-repetitive)")
    return df


def compute_genome_wide_coverage(bam_path, windows_df, threads=8, tmp_dir=None):
    """
    Compute mean coverage per 1-kb window for both unique (MAPQ≥20)
    and permissive (all primary) read sets.

    Creates a temporary filtered BAM for the unique-coverage pass then
    cleans it up.  Both passes use bedtools coverage -mean on 1-kb windows.

    Returns two numpy arrays aligned with windows_df.
    """
    import shutil

    if tmp_dir is None:
        tmp_dir = os.path.dirname(bam_path)

    # Write windows to temp BED
    tmp_bed = tempfile.NamedTemporaryFile(suffix='.bed', mode='w', delete=False, dir=tmp_dir)
    for _, row in windows_df.iterrows():
        tmp_bed.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\n")
    tmp_bed.close()

    # Also create a genome file for bedtools -sorted flag.
    # Must include EVERY contig in the BAM, in BAM sort order, or bedtools -sorted fails.
    # Derive from the BAM header (@SQ lines) instead of from AUTO_CHROMS/windows_df,
    # since the BAM may contain additional contigs (chrX, chrY, unplaced scaffolds).
    tmp_genome = tempfile.NamedTemporaryFile(suffix='.genome', mode='w', delete=False, dir=tmp_dir)
    header = subprocess.run(['samtools', 'view', '-H', bam_path],
                            capture_output=True, text=True, check=True).stdout
    for line in header.splitlines():
        if line.startswith('@SQ'):
            fields = line.split('\t')
            sn = fields[1].replace('SN:', '')
            ln = fields[2].replace('LN:', '')
            tmp_genome.write(f"{sn}\t{ln}\n")
    tmp_genome.close()

    def _run_coverage(bam_to_use, label, is_sorted_bam=True):
        """Run bedtools coverage -mean on the windows against a BAM file."""
        tmp_out = tempfile.NamedTemporaryFile(suffix='.tsv', mode='w', delete=False, dir=tmp_dir)
        tmp_out.close()

        cmd = ['bedtools', 'coverage', '-mean',
               '-a', tmp_bed.name,
               '-b', bam_to_use]
        if is_sorted_bam:
            cmd += ['-sorted', '-g', tmp_genome.name]

        print(f"  Computing {label} coverage (this may take a while)...")
        t0 = time.time()
        with open(tmp_out.name, 'w') as out_f:
            subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, check=True)
        elapsed = time.time() - t0
        print(f"    {label}: {elapsed:.0f}s")

        cov_data = pd.read_csv(tmp_out.name, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'mean_depth'])
        os.unlink(tmp_out.name)

        # Verify alignment with windows
        if len(cov_data) != len(windows_df):
            print(f"    WARNING: expected {len(windows_df)} rows, got {len(cov_data)}")

        return cov_data['mean_depth'].values

    # --- Permissive coverage (use existing BAM, which was created with -F 0x900) ---
    permissive_depths = _run_coverage(bam_path, 'permissive (all primary)')

    # --- Unique coverage: create temporary MAPQ≥20 filtered BAM ---
    tmp_unique_bam = os.path.join(tmp_dir,
                                   f".tmp_unique_mapq20_{os.path.basename(bam_path)}")
    print(f"  Creating MAPQ≥20 filtered BAM: {tmp_unique_bam} ...")
    t0 = time.time()

    # samtools view -F 0x900: keep only primary; -q 20: MAPQ >= 20
    subprocess.run(['samtools', 'view', '-@', str(max(1, threads // 2)),
                     '-F', '0x900', '-q', '20', '-b',
                     '-o', tmp_unique_bam, bam_path], check=True)
    subprocess.run(['samtools', 'index', '-@', str(max(1, threads // 2)),
                     tmp_unique_bam], check=True)
    print(f"    Filtered BAM ready in {time.time() - t0:.0f}s")

    unique_depths = _run_coverage(tmp_unique_bam, 'unique (MAPQ≥20)')

    # --- Clean up ---
    print("  Cleaning up temporary files...")
    for f in [tmp_unique_bam, tmp_unique_bam + '.bai',
              tmp_bed.name, tmp_genome.name]:
        try:
            os.unlink(f)
        except OSError:
            pass

    return unique_depths, permissive_depths


def estimate_modal_coverage(unique_depths, permissive_depths, windows_df):
    """
    Estimate genome-wide modal coverage C from non-repetitive windows.
    C is estimated from permissive coverage on repeat-free autosomal windows.
    """
    mask = ~windows_df['overlaps_repeat'].values
    mask &= np.isfinite(permissive_depths)
    mask &= (permissive_depths > 0)
    mask &= np.isfinite(unique_depths)

    dep_perm = permissive_depths[mask]
    dep_uniq = unique_depths[mask]

    C_permissive = find_coverage_mode(dep_perm)
    C_unique = find_coverage_mode(dep_uniq)

    print(f"\n  Modal permissive coverage: {C_permissive:.2f}x")
    print(f"  Modal unique (MAPQ≥20) coverage: {C_unique:.2f}x")
    print(f"  Non-repetitive windows used: {mask.sum():,}")

    return C_permissive, C_unique


# ============================================================
# Build gene regions from paralog families
# ============================================================
def build_gene_regions(paralog_df, chrom_lengths, flank_bp=25000):
    """
    From paralog_families.tsv, build region definitions for each gene.
    Returns a DataFrame with columns:
        family, gene_name, gene_id, gene_type, chrom, strand,
        gene_start, gene_end, gene_length,
        upstream_start, upstream_end,
        downstream_start, downstream_end,
        combined_start, combined_end,
        paralog_type, copy_num
    """
    regions = []
    for _, row in paralog_df.iterrows():
        chrom = row['chrom']
        chr_len = chrom_lengths.get(chrom)
        if chr_len is None:
            continue  # skip genes on chromosomes not in the assembly

        gene_start = int(row['start'])
        gene_end = int(row['end'])

        # Flanks: 25 kb up/downstream, clamped to chromosome bounds
        if row['strand'] == '+':
            upstream_start = max(0, gene_start - flank_bp)
            upstream_end = gene_start
            downstream_start = gene_end
            downstream_end = min(chr_len, gene_end + flank_bp)
        else:
            # For minus-strand genes, "upstream" is the 3' side (higher coords)
            upstream_start = gene_end
            upstream_end = min(chr_len, gene_end + flank_bp)
            downstream_start = max(0, gene_start - flank_bp)
            downstream_end = gene_start

        combined_start = max(0, gene_start - flank_bp)
        combined_end = min(chr_len, gene_end + flank_bp)

        regions.append({
            'family': row['family'],
            'gene_name': row['gene_name'],
            'gene_id': row['gene_id'],
            'gene_type': row['gene_type'],
            'chrom': chrom,
            'strand': row['strand'],
            'gene_start': gene_start,
            'gene_end': gene_end,
            'gene_length': abs(gene_end - gene_start),
            'upstream_start': upstream_start,
            'upstream_end': upstream_end,
            'upstream_length': abs(upstream_end - upstream_start),
            'downstream_start': downstream_start,
            'downstream_end': downstream_end,
            'downstream_length': abs(downstream_end - downstream_start),
            'combined_start': combined_start,
            'combined_end': combined_end,
            'combined_length': combined_end - combined_start,
            'paralog_type': row.get('paralog_type', ''),
            'copy_num': row.get('copy_num', 0),
        })

    df = pd.DataFrame(regions)
    n_auto = df[df['chrom'].isin(AUTO_CHROMS)].shape[0]
    n_sex = df[df['chrom'].isin(('chrX', 'chrY'))].shape[0]
    n_other = len(df) - n_auto - n_sex
    print(f"  Built {len(df)} gene regions from {df['family'].nunique()} families")
    print(f"    Autosomes: {n_auto}  |  Sex chr (X/Y): {n_sex}  |  Other: {n_other}")
    print(f"    Parents: {(df['gene_type'] == 'parent').sum()}")
    print(f"    Paralogs: {(df['gene_type'] != 'parent').sum()}")
    return df


# ============================================================
# Per-base coverage extraction
# ============================================================
def get_per_base_coverage(bam, chrom, start, end, min_mapq=0):
    """
    Extract per-base coverage for a region using pysam pileup.
    Returns a dict: {pos: depth} for positions with depth > 0.
    Positions with 0 coverage are NOT in the dict.
    min_mapq=0 means all reads; min_mapq>=20 filters for uniquely mapping reads.
    """
    coverage = {}
    try:
        for pileup_col in bam.pileup(
            chrom, start, end,
            truncate=True,
            min_mapping_quality=min_mapq,
            ignore_overlaps=False,
            flag_filter=0x900,  # exclude supplementary + secondary
            max_depth=5000,
        ):
            pos = pileup_col.reference_pos
            depth = pileup_col.nsegments
            if depth > 0:
                coverage[pos] = depth
    except (ValueError, OSError):
        # Region not in BAM index — return empty
        pass
    return coverage


def region_stats(cov_dict, start, end):
    """
    Compute statistics from a coverage dict {pos: depth} for region [start, end).
    Returns: median_depth, mean_depth, pct_zero, n_positions
    """
    region_len = end - start
    if region_len <= 0:
        return np.nan, np.nan, np.nan, 0

    depths = []
    zero_count = 0
    for pos in range(start, end):
        d = cov_dict.get(pos, 0)
        if d == 0:
            zero_count += 1
        else:
            depths.append(d)

    if not depths:
        return 0.0, 0.0, 100.0, region_len

    return (
        float(np.median(depths)),
        float(np.mean(depths)),
        zero_count / region_len * 100,
        region_len,
    )


def analyze_gene_coverage(bam_path, region_df, C_permissive, C_unique, threads=8):
    """
    For each gene in region_df, compute unique and permissive coverage
    across gene body, upstream, downstream, and combined regions.

    Optimised: processes chromosome by chromosome — one pileup pass per
    chromosome per MAPQ setting — then extracts per-gene statistics.
    """
    bam = pysam.AlignmentFile(bam_path, 'rb', threads=threads)

    results = []
    n_total = len(region_df)

    # Group regions by chromosome
    for chrom, group in region_df.groupby('chrom'):
        n_chr = len(group)
        print(f"  Processing {chrom}: {n_chr} genes...")

        # Determine the range covering all genes on this chromosome
        chr_c_start = int(group['combined_start'].min())
        chr_c_end = int(group['combined_end'].max())

        # --- One pileup pass per MAPQ setting ---
        # Unique (MAPQ >= 20)
        t0 = time.time()
        cov_uniq_full = get_per_base_coverage(
            bam, chrom, chr_c_start, chr_c_end, min_mapq=20)
        t_uniq = time.time() - t0

        # Permissive (all MAPQ)
        t0 = time.time()
        cov_perm_full = get_per_base_coverage(
            bam, chrom, chr_c_start, chr_c_end, min_mapq=0)
        t_perm = time.time() - t0

        print(f"    pileup: unique={t_uniq:.1f}s, permissive={t_perm:.1f}s")

        # --- Extract stats per gene ---
        for _, row in group.iterrows():
            g_start = int(row['gene_start'])
            g_end = int(row['gene_end'])
            u_start = int(row['upstream_start'])
            u_end = int(row['upstream_end'])
            d_start = int(row['downstream_start'])
            d_end = int(row['downstream_end'])
            c_start = int(row['combined_start'])
            c_end = int(row['combined_end'])

            regions_def = {
                'gene_body': (g_start, g_end),
                'upstream': (u_start, u_end),
                'downstream': (d_start, d_end),
                'combined': (c_start, c_end),
            }

            row_result = {
                'family': row['family'],
                'gene_name': row['gene_name'],
                'gene_id': row['gene_id'],
                'gene_type': row['gene_type'],
                'chrom': chrom,
                'strand': row['strand'],
                'gene_start': int(row['gene_start']),
                'gene_end': int(row['gene_end']),
                'gene_length': row['gene_length'],
                'combined_start': int(row['combined_start']),
                'combined_end': int(row['combined_end']),
                'upstream_start': int(row['upstream_start']),
                'upstream_end': int(row['upstream_end']),
                'downstream_start': int(row['downstream_start']),
                'downstream_end': int(row['downstream_end']),
                'paralog_type': row['paralog_type'],
                'copy_num': row['copy_num'],
            }

            for region_name, (r_start, r_end) in regions_def.items():
                # Unique
                med_u, mean_u, pct0_u, n_u = region_stats(
                    cov_uniq_full, r_start, r_end)
                row_result[f'{region_name}_unique_median'] = med_u
                row_result[f'{region_name}_unique_mean'] = mean_u
                row_result[f'{region_name}_unique_pct_zero'] = pct0_u
                row_result[f'{region_name}_unique_n_bases'] = n_u

                # Permissive
                med_p, mean_p, pct0_p, n_p = region_stats(
                    cov_perm_full, r_start, r_end)
                row_result[f'{region_name}_permissive_median'] = med_p
                row_result[f'{region_name}_permissive_mean'] = mean_p
                row_result[f'{region_name}_permissive_pct_zero'] = pct0_p
                row_result[f'{region_name}_permissive_n_bases'] = n_p

                # Delta: permissive - unique
                if np.isfinite(med_p) and np.isfinite(med_u):
                    row_result[f'{region_name}_median_delta'] = med_p - med_u
                else:
                    row_result[f'{region_name}_median_delta'] = np.nan

            # Normalized values (relative to modal coverage C)
            for region_name in ['gene_body', 'upstream', 'downstream', 'combined']:
                for track in ['unique', 'permissive']:
                    med_key = f'{region_name}_{track}_median'
                    norm_key = f'{region_name}_{track}_norm_median'
                    C_ref = C_unique if track == 'unique' else C_permissive
                    med_val = row_result[med_key]
                    if np.isfinite(med_val) and np.isfinite(C_ref) and C_ref > 0:
                        row_result[norm_key] = med_val / C_ref
                    else:
                        row_result[norm_key] = np.nan

            results.append(row_result)

    bam.close()
    result_df = pd.DataFrame(results)
    print(f"  Coverage computed for {len(result_df)} genes")
    return result_df


# ============================================================
# Pair-level classification
# ============================================================
def classify_pairs(cov_df):
    """
    For each family, classify the parent-paralog relationship.
    Returns a list of per-pair classification records.
    """
    pairs = []

    for family, group in cov_df.groupby('family'):
        parents = group[group['gene_type'] == 'parent']
        paralogs = group[group['gene_type'] != 'parent']

        if len(parents) == 0 or len(paralogs) == 0:
            continue

        for _, par_row in paralogs.iterrows():
            # For simplicity, use the first parent
            par_row_parent = parents.iloc[0]

            pair_record = _classify_single_pair(par_row_parent, par_row)
            pairs.append(pair_record)

    pair_df = pd.DataFrame(pairs)
    return pair_df


def _classify_single_pair(parent_row, paralog_row):
    """
    Classify a single parent–paralog pair according to the screening logic.
    """
    pair = {
        'family': parent_row['family'],
        'parent_gene': parent_row['gene_name'],
        'paralog_gene': paralog_row['gene_name'],
        'parent_chrom': parent_row['chrom'],
        'paralog_chrom': paralog_row['chrom'],
        'parent_type': parent_row['gene_type'],
        'paralog_type': paralog_row['gene_type'],
        'same_chromosome': parent_row['chrom'] == paralog_row['chrom'],
        'on_sex_chromosome': (parent_row['chrom'] in ('chrX', 'chrY')
                              or paralog_row['chrom'] in ('chrX', 'chrY')),
    }

    # Extract key metrics
    P_UF_unique = parent_row['combined_unique_norm_median']
    P_UF_perm = parent_row['combined_permissive_norm_median']
    P_FL_uq = parent_row['upstream_unique_norm_median']
    P_FL_dq = parent_row['downstream_unique_norm_median']
    P_FL_up = parent_row['upstream_permissive_norm_median']
    P_FL_dp = parent_row['downstream_permissive_norm_median']
    P_GB_uq = parent_row['gene_body_unique_norm_median']
    P_GB_perm = parent_row['gene_body_permissive_norm_median']

    C_UF_unique = paralog_row['combined_unique_norm_median']
    C_UF_perm = paralog_row['combined_permissive_norm_median']
    C_FL_uq = paralog_row['upstream_unique_norm_median']
    C_FL_dq = paralog_row['downstream_unique_norm_median']
    C_FL_up = paralog_row['upstream_permissive_norm_median']
    C_FL_dp = paralog_row['downstream_permissive_norm_median']
    C_GB_uq = paralog_row['gene_body_unique_norm_median']
    C_GB_perm = paralog_row['gene_body_permissive_norm_median']

    # Flank mappability check: fraction of zero-coverage bases in unique track
    P_FL_up_pct0 = parent_row['upstream_unique_pct_zero']
    P_FL_dn_pct0 = parent_row['downstream_unique_pct_zero']
    C_FL_up_pct0 = paralog_row['upstream_unique_pct_zero']
    C_FL_dn_pct0 = paralog_row['downstream_unique_pct_zero']

    # Combined unique coverage (for sum check)
    sum_uq_med = np.nansum([P_UF_unique, C_UF_unique])

    # Store all metrics
    for k in ['unique_median', 'unique_mean', 'unique_norm_median',
              'permissive_median', 'permissive_mean', 'permissive_norm_median',
              'unique_pct_zero', 'permissive_pct_zero',
              'median_delta']:
        for region in ['gene_body', 'upstream', 'downstream', 'combined']:
            pair[f'parent_{region}_{k}'] = parent_row.get(f'{region}_{k}', np.nan)
            pair[f'paralog_{region}_{k}'] = paralog_row.get(f'{region}_{k}', np.nan)

    # --- Classification logic ---

    # Tolerance for "near C": within 20% of C (0.8–1.2)
    # For "near C/2": within 0.3–0.7 (relaxed)
    def near_diploid(val):
        return (not np.isnan(val)) and 0.65 <= val <= 1.35

    def near_haploid(val):
        return (not np.isnan(val)) and 0.25 <= val <= 0.75

    def is_informative(pct0):
        """Flank is informative if < 50% zero-coverage bases."""
        return (not np.isnan(pct0)) and pct0 < 50

    # Summarize flank states
    parent_flank_near_C = near_diploid(P_FL_uq) and near_diploid(P_FL_dq)
    paralog_flank_near_C = near_diploid(C_FL_uq) and near_diploid(C_FL_dq)

    parent_flank_near_half = near_haploid(P_FL_uq) or near_haploid(P_FL_dq)
    paralog_flank_near_half = near_haploid(C_FL_uq) or near_haploid(C_FL_dq)

    parent_flank_info = is_informative(P_FL_up_pct0) and is_informative(P_FL_dn_pct0)
    paralog_flank_info = is_informative(C_FL_up_pct0) and is_informative(C_FL_dn_pct0)

    # Gene body states
    parent_gb_low = (not np.isnan(P_GB_uq)) and P_GB_uq < 0.65
    paralog_gb_low = (not np.isnan(C_GB_uq)) and C_GB_uq < 0.65

    parent_gb_perm_ok = (not np.isnan(P_GB_perm)) and P_GB_perm > 0.65
    paralog_gb_perm_ok = (not np.isnan(C_GB_perm)) and C_GB_perm > 0.65

    # --- Decision tree ---

    # Rule 1: Likely genuine duplication
    # Both loci have unique flank coverage near C
    if parent_flank_near_C and paralog_flank_near_C:
        pair['classification'] = CLASS_GENUINE
        pair['classification_reason'] = (
            'Both parent and paralog unique flanks near diploid coverage C'
        )

    # Rule 3: Likely recent paralog with mapping ambiguity
    # Gene-body unique coverage is low but permissive is OK, flanks at C
    elif ((parent_gb_low and parent_gb_perm_ok) or (paralog_gb_low and paralog_gb_perm_ok)) \
            and parent_flank_near_C and paralog_flank_near_C:
        pair['classification'] = CLASS_AMBIG
        pair['classification_reason'] = (
            'Gene body unique coverage low with OK permissive; flanks at C'
        )

    # Rule 2: Suspicious for haplotypic duplication
    # Both loci at ~C/2, including flank evidence, and sums to ~C
    elif (near_haploid(P_UF_unique) and near_haploid(C_UF_unique)
          and near_diploid(sum_uq_med)
          and parent_flank_near_half and paralog_flank_near_half
          and parent_flank_info and paralog_flank_info):
        pair['classification'] = CLASS_HAPLOTIG
        pair['classification_reason'] = (
            'Both loci near C/2 with half-depth flanks; combined coverage near C'
        )

    # Softer haplotig call: half-depth across gene + flanks but sum check borderline
    elif (near_haploid(P_UF_unique) and near_haploid(C_UF_unique)
          and parent_flank_near_half and paralog_flank_near_half
          and (parent_flank_info or paralog_flank_info)):
        pair['classification'] = CLASS_HAPLOTIG
        pair['classification_reason'] = (
            'Both loci near C/2 with half-depth signal in at least one flank pair'
        )

    # Rule 4: Inconclusive
    else:
        pair['classification'] = CLASS_INCONCLUSIVE
        reasons = []
        if not parent_flank_info and not paralog_flank_info:
            reasons.append('insufficient mappable flank sequence at both loci')
        if not reasons:
            reasons.append('coverage pattern does not clearly match any category')
        pair['classification_reason'] = '; '.join(reasons)

    # Ploidy caveat: chrX/Y are haploid in males, so C/2 is the expected
    # coverage for sex-chromosome genes, not necessarily a haplotig signal.
    if pair['on_sex_chromosome']:
        pair['classification_reason'] += ' [NOTE: sex chromosome — haploid baseline expected]'

    return pair


# ============================================================
# Plot generation
# ============================================================
def make_coverage_plot(parent_row, paralog_row, pair_classification,
                        cov_df, bam, out_dir, flank_bp, C_permissive, C_unique):
    """
    Generate a normalized coverage plot for a single parent–paralog pair
    across gene body ± flank region.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Open BAM handle (caller manages open/close).
    """
    fig, axes = plt.subplots(1, 2, figsize=(20, 7))

    for ax_idx, (gene_row, label) in enumerate([
        (parent_row, 'Parent'),
        (paralog_row, 'Paralog'),
    ]):
        ax = axes[ax_idx]

        chrom = gene_row['chrom']
        g_start = int(gene_row['gene_start'])
        g_end = int(gene_row['gene_end'])
        c_start = int(gene_row['combined_start'])
        c_end = int(gene_row['combined_end'])
        strand = gene_row['strand']

        # Get per-base coverage
        cov_uniq = get_per_base_coverage(bam, chrom, c_start, c_end, min_mapq=20)
        cov_perm = get_per_base_coverage(bam, chrom, c_start, c_end, min_mapq=0)

        # Build position arrays with normalized values
        positions_rel = np.arange(c_start, c_end) - g_start  # relative to gene start
        x_vals = positions_rel / 1000  # kb

        uniq_vals = np.array([cov_uniq.get(p, 0) for p in range(c_start, c_end)],
                              dtype=float)
        perm_vals = np.array([cov_perm.get(p, 0) for p in range(c_start, c_end)],
                              dtype=float)

        # Normalize by modal coverage
        uniq_norm = uniq_vals / C_unique if C_unique > 0 else uniq_vals
        perm_norm = perm_vals / C_permissive if C_permissive > 0 else perm_vals

        # Plot permissive (lighter)
        ax.fill_between(x_vals, 0, perm_norm, alpha=0.15,
                         color='steelblue', label='Permissive (all MAPQ)')
        ax.plot(x_vals, perm_norm, color='steelblue', linewidth=0.5, alpha=0.7)

        # Plot unique (darker)
        ax.fill_between(x_vals, 0, uniq_norm, alpha=0.2,
                         color='darkorange', label='Unique (MAPQ≥20)')
        ax.plot(x_vals, uniq_norm, color='darkorange', linewidth=0.5, alpha=0.8)

        # Gene body highlight
        gene_x_start = 0
        gene_x_end = (g_end - g_start) / 1000
        ax.axvspan(gene_x_start, gene_x_end, alpha=0.08, color='red')
        ax.axvline(x=gene_x_start, color='red', linestyle='--', linewidth=0.8, alpha=0.6)
        ax.axvline(x=gene_x_end, color='red', linestyle='--', linewidth=0.8, alpha=0.6)

        # Diploid reference line
        ax.axhline(y=1.0, color='green', linestyle='--', linewidth=1.0, alpha=0.5,
                    label='C (diploid)')
        ax.axhline(y=0.5, color='red', linestyle=':', linewidth=1.0, alpha=0.4,
                    label='C/2 (haploid)')

        ax.set_xlabel('Position relative to gene start (kb)', fontsize=11)
        ax.set_ylabel(f'Normalized coverage (× / C)', fontsize=11)
        gene_id = gene_row.get('gene_id', gene_row['gene_name'])
        ax.set_title(f'{label}: {gene_row["gene_name"]}\n'
                      f'{chrom}:{g_start:,}-{g_end:,}  ({strand})',
                      fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='upper right', ncol=2)

        # Set y-limit
        max_y = max(np.percentile(perm_norm[perm_norm > 0], 99) if np.any(perm_norm > 0) else 2, 2.5)
        ax.set_ylim(0, min(max_y * 1.2, 5))
        ax.set_xlim(x_vals[0], x_vals[-1])

        # Annotate median values
        for region_name, r_start, r_end in [
            ('upstream', int(gene_row['upstream_start']), int(gene_row['upstream_end'])),
            ('gene', g_start, g_end),
            ('downstream', int(gene_row['downstream_start']), int(gene_row['downstream_end'])),
        ]:
            r_uniq = np.array([cov_uniq.get(p, 0) for p in range(r_start, r_end)], dtype=float)
            r_perm = np.array([cov_perm.get(p, 0) for p in range(r_start, r_end)], dtype=float)
            med_u = np.median(r_uniq) / C_unique if C_unique and len(r_uniq) > 0 else np.nan
            med_p = np.median(r_perm) / C_permissive if C_permissive and len(r_perm) > 0 else np.nan

            if not np.isnan(med_u):
                mid_x = (r_start + r_end) / 2 - g_start
                ax.annotate(f'u:{med_u:.2f}\np:{med_p:.2f}',
                            xy=(mid_x / 1000, 0.95 * ax.get_ylim()[1]),
                            fontsize=7, ha='center', va='top',
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

    # Overall title
    fig.suptitle(f'{parent_row["family"]}: {parent_row["gene_name"]} — {paralog_row["gene_name"]}\n'
                  f'Classification: {pair_classification}',
                  fontsize=14, fontweight='bold', y=1.01)

    plt.tight_layout()
    safe_family = parent_row['family'].replace('/', '_').replace(' ', '_')
    safe_parent = parent_row['gene_name'].replace('/', '_')
    safe_paralog = paralog_row['gene_name'].replace('/', '_')
    plot_path = os.path.join(out_dir, 'plots',
                              f'{safe_family}__{safe_parent}__{safe_paralog}.png')
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    return plot_path


def make_summary_figure(pair_df, out_dir):
    """Make overview summary figure of classification results."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Panel A: Classification counts
    ax = axes[0, 0]
    class_counts = pair_df['classification'].value_counts()
    class_order = [CLASS_GENUINE, CLASS_HAPLOTIG, CLASS_AMBIG, CLASS_INCONCLUSIVE]
    colors = {'likely_genuine_duplication': '#2ca02c',
              'suspicious_haplotig': '#d62728',
              'recent_paralog_mapping_ambiguity': '#ff7f0e',
              'inconclusive': '#7f7f7f'}
    present = [c for c in class_order if c in class_counts.index]
    counts = [class_counts.get(c, 0) for c in present]
    bar_colors = [colors.get(c, 'gray') for c in present]
    bars = ax.bar(range(len(present)), counts, color=bar_colors, edgecolor='white')
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(count), ha='center', fontsize=10, fontweight='bold')
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([c.replace('_', '\n') for c in present], fontsize=9)
    ax.set_ylabel('Number of pairs')
    ax.set_title('A. Classification Summary', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel B: Unique vs permissive combined coverage scatter by classification
    ax = axes[0, 1]
    for cls in present:
        subset = pair_df[pair_df['classification'] == cls]
        ax.scatter(subset['parent_combined_unique_norm_median'],
                   subset['paralog_combined_unique_norm_median'],
                   c=colors.get(cls, 'gray'), alpha=0.5, s=15,
                   label=f'{cls} ({len(subset)})')
    ax.axhline(y=1.0, color='green', linestyle='--', alpha=0.4)
    ax.axvline(x=1.0, color='green', linestyle='--', alpha=0.4)
    ax.axhline(y=0.5, color='red', linestyle=':', alpha=0.3)
    ax.axvline(x=0.5, color='red', linestyle=':', alpha=0.3)
    ax.set_xlabel('Parent combined unique norm. coverage')
    ax.set_ylabel('Paralog combined unique norm. coverage')
    ax.set_title('B. Parent vs Paralog Unique Coverage', fontsize=13, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel C: Upstream flank coverage distribution by class
    ax = axes[1, 0]
    flank_data = []
    flank_labels = []
    for cls in present:
        subset = pair_df[pair_df['classification'] == cls]
        vals = pd.concat([
            subset['parent_upstream_unique_norm_median'].dropna(),
            subset['paralog_upstream_unique_norm_median'].dropna(),
            subset['parent_downstream_unique_norm_median'].dropna(),
            subset['paralog_downstream_unique_norm_median'].dropna(),
        ])
        if len(vals) > 0:
            flank_data.append(vals.clip(0, 2.5).values)
            flank_labels.append(f'{cls}\n(n={len(subset)})')

    if flank_data:
        bp = ax.boxplot(flank_data, labels=flank_labels, patch_artist=True,
                         showfliers=False)
        for patch, cls in zip(bp['boxes'], present):
            patch.set_facecolor(colors.get(cls, 'gray'))
            patch.set_alpha(0.5)
        ax.axhline(y=1.0, color='green', linestyle='--', alpha=0.4)
        ax.axhline(y=0.5, color='red', linestyle=':', alpha=0.3)
        ax.set_ylabel('Unique flank norm. coverage')
        ax.set_title('C. Flank Unique Coverage by Class', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Panel D: Coverage delta (permissive - unique) by class
    ax = axes[1, 1]
    delta_data = []
    delta_labels = []
    for cls in present:
        subset = pair_df[pair_df['classification'] == cls]
        vals = subset['parent_combined_median_delta'].dropna()
        if len(vals) > 0:
            delta_data.append(vals.clip(-0.5, 5).values)
            delta_labels.append(f'{cls}\n(n={len(subset)})')

    if delta_data:
        bp2 = ax.boxplot(delta_data, labels=delta_labels, patch_artist=True,
                          showfliers=False)
        for patch, cls in zip(bp2['boxes'], present):
            patch.set_facecolor(colors.get(cls, 'gray'))
            patch.set_alpha(0.5)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        ax.set_ylabel('Coverage delta (permissive − unique)')
        ax.set_title('D. Permissive−Unique Coverage Delta', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    summary_path = os.path.join(out_dir, 'classification_summary.png')
    plt.savefig(summary_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Summary figure: {summary_path}")


# ============================================================
# Paralog families loader
# ============================================================
def load_paralog_families(tsv_path):
    """Load paralog_families.tsv and return a standardized DataFrame."""
    df = pd.read_csv(tsv_path, sep='\t')
    # Expected columns: family, gene_name, gene_id, gene_type, chrom, start, end, strand, length, paralog_type, copy_num
    required = ['family', 'gene_name', 'gene_type', 'chrom', 'start', 'end', 'strand']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column '{col}' in {tsv_path}. "
                             f"Available columns: {list(df.columns)}")
    print(f"  Loaded {len(df)} gene entries from {len(df['family'].unique())} families")
    return df


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description='Read-depth haplotig screening for parent–paralog pairs')
    parser.add_argument('--bam', required=True, help='Sorted, indexed BAM file')
    parser.add_argument('--fasta', required=True, help='Reference assembly FASTA')
    parser.add_argument('--paralog-families', required=True,
                        help='paralog_families.tsv with parent/paralog gene coordinates')
    parser.add_argument('--genome-lengths', required=True,
                        help='Two-column file: chrom\\tlength')
    parser.add_argument('--out-dir', required=True, help='Output directory')
    parser.add_argument('--repeat-bed', default=None,
                        help='Repeat annotation BED/GFF (optional)')
    parser.add_argument('--flank-bp', type=int, default=25000,
                        help='Flank size in bp (default: 25000)')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads (default: 8)')
    parser.add_argument('--max-families', type=int, default=None,
                        help='Limit to first N families for testing (default: all)')
    parser.add_argument('--skip-genome-wide', action='store_true',
                        help='Skip genome-wide coverage estimation (use provided --C)')
    parser.add_argument('--C-permissive', type=float, default=None,
                        help='Pre-computed modal permissive coverage')
    parser.add_argument('--C-unique', type=float, default=None,
                        help='Pre-computed modal unique coverage')
    parser.add_argument('--skip-plots', action='store_true',
                        help='Skip per-pair coverage plots')
    parser.add_argument('--no-cache', action='store_true',
                        help='Ignore cached output files and re-run all phases')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 70)
    print("READ-DEPTH HAPLOTIG SCREENING PIPELINE")
    print("=" * 70)
    print(f"BAM: {args.bam}")
    print(f"Assembly: {args.fasta}")
    print(f"Paralog families: {args.paralog_families}")
    print(f"Output: {args.out_dir}")
    print(f"Flank: ±{args.flank_bp} bp")
    print(f"Threads: {args.threads}")
    print()

    # ============================================================
    # Phase 0: Load inputs
    # ============================================================
    print("PHASE 0: Loading inputs...")

    # Chromosome lengths
    chr_lengths = {}
    with open(args.genome_lengths) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                chr_lengths[parts[0]] = int(parts[1])
    print(f"  Loaded {len(chr_lengths)} chromosome/scaffold lengths")
    print(f"  Autosomes (for coverage baseline): {sum(1 for c in AUTO_CHROMS if c in chr_lengths)}")

    # Paralog families
    paralog_df = load_paralog_families(args.paralog_families)
    if args.max_families:
        families = paralog_df['family'].unique()[:args.max_families]
        paralog_df = paralog_df[paralog_df['family'].isin(families)]
        print(f"  Limited to first {args.max_families} families: {len(paralog_df)} genes")

    # Build gene regions
    region_df = build_gene_regions(paralog_df, chr_lengths, args.flank_bp)

    # ============================================================
    # Phase 1: Genome-wide modal coverage
    # ============================================================
    if args.skip_genome_wide and args.C_permissive is not None and args.C_unique is not None:
        C_permissive, C_unique = args.C_permissive, args.C_unique
        print(f"\nPHASE 1: Using provided modal coverage: "
              f"permissive={C_permissive:.2f}x, unique={C_unique:.2f}x")
    else:
        genome_wide_path = os.path.join(args.out_dir, 'genome_wide_1kb_coverage.tsv')
        if os.path.exists(genome_wide_path) and not args.no_cache:
            print("\nPHASE 1: Loading cached genome-wide coverage...")
            windows_df = pd.read_csv(genome_wide_path, sep='\t')
            unique_depths = windows_df['cov_unique'].values
            permissive_depths = windows_df['cov_permissive'].values
            windows_df = windows_df.drop(columns=['cov_unique', 'cov_permissive'])
            C_permissive, C_unique = estimate_modal_coverage(
                unique_depths, permissive_depths, windows_df)
            print(f"  Loaded {len(windows_df):,} windows from cache")
        else:
            print("\nPHASE 1: Estimating genome-wide modal coverage...")

            # Load repeat mask (if provided)
            repeats = {}
            if args.repeat_bed and os.path.exists(args.repeat_bed):
                print(f"  Loading repeat annotation: {args.repeat_bed}")
                repeats = load_repeat_mask(args.repeat_bed)
                n_intervals = sum(len(v) for v in repeats.values())
                print(f"  Loaded {n_intervals:,} repeat intervals across {len(repeats)} contigs")

            # Build windows
            windows_df = build_windows(chr_lengths, repeats if repeats else None)

            # Compute coverage per window
            t0 = time.time()
            unique_depths, permissive_depths = compute_genome_wide_coverage(
                args.bam, windows_df, args.threads)
            print(f"  Window coverage done in {time.time() - t0:.0f}s")

            # Estimate modal coverage
            C_permissive, C_unique = estimate_modal_coverage(
                unique_depths, permissive_depths, windows_df)

            # Save window-level coverage for reference
            windows_df['cov_unique'] = unique_depths
            windows_df['cov_permissive'] = permissive_depths
            windows_df.to_csv(genome_wide_path, sep='\t', index=False)
            print(f"  Saved genome-wide window coverage")

    if np.isnan(C_permissive) or np.isnan(C_unique):
        print("ERROR: Could not estimate modal coverage. Check BAM and genome length files.")
        sys.exit(1)

    print(f"\n  Modal permissive coverage C = {C_permissive:.2f}x")
    print(f"  Modal unique coverage C = {C_unique:.2f}x")

    # Also plot genome-wide coverage histogram
    if not args.skip_genome_wide:
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            mask = np.isfinite(permissive_depths) & (permissive_depths > 0)
            dep_plot = permissive_depths[mask]
            if len(dep_plot) > 0:
                lo, hi = np.percentile(dep_plot, [0.1, 99.9])
                ax.hist(dep_plot.clip(lo, hi), bins=100, color='steelblue',
                         alpha=0.7, edgecolor='white', density=True)
                ax.axvline(x=C_permissive, color='red', linestyle='--', linewidth=2,
                            label=f'Mode = {C_permissive:.2f}x')
                ax.set_xlabel('Mean coverage per 1-kb window', fontsize=12)
                ax.set_ylabel('Density', fontsize=12)
                ax.set_title('Genome-wide Coverage Distribution (Permissive)', fontsize=14)
                ax.legend()
                fig.savefig(os.path.join(args.out_dir, 'genome_wide_coverage_hist.png'),
                            dpi=150, bbox_inches='tight')
                plt.close()
        except Exception as e:
            print(f"  Warning: Could not generate coverage histogram: {e}")

    # ============================================================
    # Phase 2: Per-gene coverage analysis
    # ============================================================
    cov_path = os.path.join(args.out_dir, 'per_gene_coverage.tsv')
    if os.path.exists(cov_path) and not args.no_cache:
        print("\nPHASE 2: Loading cached per-gene coverage...")
        cov_df = pd.read_csv(cov_path, sep='\t')
        print(f"  Loaded {len(cov_df)} gene entries from cache")
    else:
        print("\nPHASE 2: Computing per-gene coverage...")
        t0 = time.time()

        cov_df = analyze_gene_coverage(
            args.bam, region_df, C_permissive, C_unique, args.threads)

        print(f"  Done in {time.time() - t0:.0f}s")

        # Save per-gene coverage
        cov_df.to_csv(cov_path, sep='\t', index=False)
        print(f"  Saved: {cov_path}")

    # ============================================================
    # Phase 3: Pair classification
    # ============================================================
    pair_path = os.path.join(args.out_dir, 'pair_classification.tsv')
    suspicious_path = os.path.join(args.out_dir, 'suspicious_pairs_manual_review.tsv')
    if os.path.exists(pair_path) and not args.no_cache:
        print("\nPHASE 3: Loading cached pair classification...")
        pair_df = pd.read_csv(pair_path, sep='\t')
        suspicious = pair_df[pair_df['classification'] == CLASS_HAPLOTIG].copy()
        print(f"  Loaded {len(pair_df)} pairs from cache")
    else:
        print("\nPHASE 3: Classifying parent–paralog pairs...")
        pair_df = classify_pairs(cov_df)

        # Save pair classification
        pair_df.to_csv(pair_path, sep='\t', index=False)
        print(f"\n  Saved: {pair_path}")

        # Suspicious pairs for manual review
        suspicious = pair_df[pair_df['classification'] == CLASS_HAPLOTIG].copy()
        suspicious.to_csv(suspicious_path, sep='\t', index=False)

    # Summary
    class_counts = pair_df['classification'].value_counts()
    print(f"\n  Classification results:")
    for cls, count in class_counts.items():
        pct = count / len(pair_df) * 100
        print(f"    {cls}: {count} ({pct:.1f}%)")

    print(f"  Suspicious pairs requiring manual review: {len(suspicious)}")
    if len(suspicious) > 0:
        print(f"\n  Suspicious pairs:")
        for _, row in suspicious.iterrows():
            print(f"    {row['family']}: {row['parent_gene']} — {row['paralog_gene']}")
            print(f"      {row['classification_reason']}")

    # ============================================================
    # Phase 4: Plot generation
    # ============================================================
    summary_fig_path = os.path.join(args.out_dir, 'classification_summary.png')
    if os.path.exists(summary_fig_path) and not args.no_cache:
        print("\nPHASE 4: Cached summary figure exists, skipping plot generation.")
        print(f"  Use --no-cache to force regeneration.")
    else:
        print("\nPHASE 4: Generating plots...")

        # Summary figure
        make_summary_figure(pair_df, args.out_dir)

    # Per-pair plots (for suspicious + a sampling of others)
    if not args.skip_plots:
        plots_dir = os.path.join(args.out_dir, 'plots')
        existing_plots = (os.path.isdir(plots_dir) and
                          len([f for f in os.listdir(plots_dir) if f.endswith('.png')]) > 0)
        if existing_plots and not args.no_cache:
            n_existing = len([f for f in os.listdir(plots_dir) if f.endswith('.png')])
            print(f"  Found {n_existing} cached per-pair plots, skipping.")
        else:
            # Always plot suspicious pairs
            plot_pairs = pair_df[pair_df['classification'] == CLASS_HAPLOTIG]

            # Add a sample of other classes (up to 50 each)
            for cls in [CLASS_GENUINE, CLASS_AMBIG, CLASS_INCONCLUSIVE]:
                cls_subset = pair_df[pair_df['classification'] == cls]
                if len(cls_subset) > 0:
                    sample_n = min(50, len(cls_subset))
                    plot_pairs = pd.concat([plot_pairs, cls_subset.sample(
                        sample_n, random_state=42)])

            print(f"  Generating {len(plot_pairs)} per-pair coverage plots...")
            n_plotted = 0
            bam_plot = pysam.AlignmentFile(args.bam, 'rb')

            for _, pair_row in plot_pairs.iterrows():
                # Get parent and paralog rows from cov_df
                parent_mask = (cov_df['gene_name'] == pair_row['parent_gene'])
                paralog_mask = (cov_df['gene_name'] == pair_row['paralog_gene'])

                parent_row = cov_df[parent_mask].iloc[0] if parent_mask.any() else None
                paralog_row = cov_df[paralog_mask].iloc[0] if paralog_mask.any() else None

                if parent_row is None or paralog_row is None:
                    continue

                try:
                    make_coverage_plot(
                        parent_row, paralog_row,
                        pair_row['classification'],
                        cov_df, bam_plot, args.out_dir,
                        args.flank_bp, C_permissive, C_unique)
                    n_plotted += 1
                except Exception as e:
                    print(f"    Warning: Could not plot {pair_row['family']}: {e}")

            bam_plot.close()
            print(f"  Generated {n_plotted} plots")

    # ============================================================
    # Done
    # ============================================================
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"\nOutput directory: {args.out_dir}")
    print(f"Files:")
    for f in sorted(os.listdir(args.out_dir)):
        fpath = os.path.join(args.out_dir, f)
        if os.path.isfile(fpath):
            print(f"  {f}")
    plots_dir = os.path.join(args.out_dir, 'plots')
    if os.path.isdir(plots_dir):
        for f in sorted(os.listdir(plots_dir)):
            if f.endswith('.png'):
                print(f"  plots/{f}")
                break
        print(f"  plots/ (total: {len(os.listdir(plots_dir))})")


if __name__ == '__main__':
    main()
