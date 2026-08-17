#!/usr/bin/env python3
"""Build the v4 BUSCO-calibrated locus-centric classifier category exploration notebook.

Panel structure (each saved as separate file for Illustrator assembly):
  A. Locus-level coverage classification criteria (ax.table)
  B. Distribution of locus coverage categories (pie chart)
  C. Coverage category by paralog type (stacked bar %)
  D-I. Representative locus examples — one paralog per category
  J. Pair-level interpretation framework (ax.table)

V4 uses thresholds calibrated from 12,057 complete autosomal BUSCO loci
(central 95% reference interval [0.636, 1.394]).
"""

import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {
        "display_name": "Python (genome-assembly)",
        "language": "python",
        "name": "genome-assembly",
    },
    "language_info": {
        "name": "python",
        "version": "3.10.0",
    },
}

cells = []


def add_md(source):
    cells.append(nbf.v4.new_markdown_cell(source))


def add_code(source):
    cells.append(nbf.v4.new_code_cell(source))


# ===========================================================================
# Title
# ===========================================================================
add_md("""# BUSCO-calibrated HiFi read-depth evaluation of annotated paralog loci

**Purpose**: This notebook evaluates HiFi read depth first at the level of individual
paralog loci (Stage 1), then integrates locus-level observations into pair-level
biological interpretations (Stage 2).

**Data sources**:
- `read-depth-screen-v4-busco-calibrated/locus_classification_v4_busco.tsv` — 3,347 loci
- `read-depth-screen-v4-busco-calibrated/pair_integration_v4_busco.tsv` — 2,170 pairs
- `paralog_genomic_summary.tsv` — paralog metadata (type, expression, genomic context)

**V4 thresholds** (BUSCO-calibrated from 12,057 complete autosomal single-copy loci):
- Expected depth: [0.636, 1.394] — central 95% of autosomal BUSCO
- Half-depth-like: [0.30, 0.55] — false half-depth rate 0.6% among BUSCO
- Very-low: [0.20, 0.30) — below empirical BUSCO range (min = 0.394)
- Intermediate: (0.55, 0.636)
""")

# ===========================================================================
# Section 1: Load and prepare data
# ===========================================================================
add_md("## 1. Load and prepare data")

add_code("""import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')
import re
import textwrap

# ---------------------------------------------------------------------------
# Publication-style settings
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'figure.dpi': 150,
    'font.family': 'sans-serif',
    'font.size': 9,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'legend.fontsize': 8,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# ---------------------------------------------------------------------------
# Color palettes
# ---------------------------------------------------------------------------
LOCUS_COLORS = {
    'normal_depth': '#2ECC40',
    'consistent_half_depth': '#E74C3C',
    'very_low_coverage': '#FF4136',
    'localized_or_asymmetric_low_depth': '#FF851B',
    'unique_mapping_deficit': '#9B59B6',
    'inconclusive': '#AAAAAA',
}

PAIR_COLORS = {
    'both_expected_depth': '#2ECC40',
    'one_expected_one_atypical': '#FFDC00',
    'both_atypical': '#F39C12',
    'asymmetric_coverage': '#FF851B',
    'unresolved': '#AAAAAA',
    'concordant_half_depth_sum_valid': '#E74C3C',
}

PARALOG_TYPE_COLORS = {
    'l': '#3498DB',
    'dl': '#E67E22',
    'rl': '#2ECC40',
}

# ---------------------------------------------------------------------------
# Display-name mappings (publication-facing labels)
# Internal names are preserved in code and TSV outputs; these are used
# only for figures, tables, legends, and panel annotations.
# ---------------------------------------------------------------------------
LOCUS_DISPLAY = {
    'normal_depth': 'EXPECTED DEPTH',
    'consistent_half_depth': 'HALF-DEPTH',
    'very_low_coverage': 'VERY LOW COVERAGE',
    'localized_or_asymmetric_low_depth': 'LOCALIZED LOW DEPTH',
    'unique_mapping_deficit': 'MAPPING AMBIGUITY',
    'inconclusive': 'UNEVALUABLE',
}

PAIR_DISPLAY = {
    'concordant_half_depth_sum_valid': 'CONCORDANT HALF-DEPTH',
    'both_expected_depth': 'BOTH LOCI AT EXPECTED DEPTH',
    'asymmetric_coverage': 'ASYMMETRIC COVERAGE',
    'one_expected_one_atypical': 'ONE EXPECTED, ONE ATYPICAL',
    'both_atypical': 'BOTH ATYPICAL',
    'unresolved': 'INDETERMINATE',
}

LOCUS_INTERPRETATION = {
    'normal_depth': (
        'Coverage within the empirical\\n'
        'range of single-copy BUSCO\\n'
        'genes. Consistent with the\\n'
        'expected genomic representation\\n'
        'of the locus.'
    ),
    'consistent_half_depth': (
        'Persistent locus-wide half-\\n'
        'depth signal. Overlaps extreme\\n'
        'BUSCO tail (false half-depth\\n'
        'rate 0.6%). Requires pair-level\\n'
        'and genomic-context\\n'
        'interpretation.'
    ),
    'very_low_coverage': (
        'Below the empirical range of\\n'
        'all complete autosomal BUSCO\\n'
        'genes (min=0.394). More\\n'
        'consistent with mapping\\n'
        'dropout, deletions, or assembly\\n'
        'collapse than with haploid\\n'
        'copy number.'
    ),
    'localized_or_asymmetric_low_depth': (
        'Localized coverage anomaly;\\n'
        'review for mapping, CNV,\\n'
        'repeat, or assembly-boundary\\n'
        'effects.'
    ),
    'unique_mapping_deficit': (
        'Coverage is supported, but\\n'
        'locus-specific read assignment\\n'
        'is limited by high sequence\\n'
        'similarity.'
    ),
    'inconclusive': (
        'Coverage cannot be interpreted\\n'
        'confidently; manual review may\\n'
        'be required.'
    ),
}

PAIR_INTERPRETATION_TEXT = {
    'concordant_half_depth_sum_valid': (
        'Both loci exhibit concordant\\n'
        'half-depth coverage with a\\n'
        'valid sum check. Pattern is\\n'
        'consistent with a retained\\n'
        'haplotypic duplication.'
    ),
    'both_expected_depth': (
        'Both loci exhibit expected\\n'
        'diploid coverage. Pattern is\\n'
        'consistent with independently\\n'
        'assembled paralogous loci.'
    ),
    'asymmetric_coverage': (
        'The two loci exhibit different\\n'
        'coverage states. Pattern may\\n'
        'reflect copy-number differences,\\n'
        'assembly effects, or biological\\n'
        'divergence.'
    ),
    'one_expected_one_atypical': (
        'One locus exhibits expected\\n'
        'coverage while the other shows\\n'
        'an atypical coverage pattern.\\n'
        'Additional genomic evidence is\\n'
        'required.'
    ),
    'both_atypical': (
        'Both loci exhibit atypical\\n'
        'coverage patterns. Relationship\\n'
        'cannot be interpreted from\\n'
        'read depth alone.'
    ),
    'unresolved': (
        'Neither locus provides evaluable\\n'
        'coverage data (both flanks\\n'
        'uninformative or no pattern\\n'
        'matched). Coverage cannot\\n'
        'contribute to interpretation.'
    ),
}

# ---------------------------------------------------------------------------
# BUSCO-calibrated thresholds
# ---------------------------------------------------------------------------
EXPECTED_LO, EXPECTED_HI = 0.636, 1.394
HALF_LO, HALF_HI = 0.30, 0.55
VERY_LOW_LO, VERY_LOW_HI = 0.20, 0.30

# Locus category display order
LOCUS_ORDER = ['normal_depth', 'consistent_half_depth', 'very_low_coverage',
               'localized_or_asymmetric_low_depth', 'unique_mapping_deficit',
               'inconclusive']

# Pair interpretation order (by count, most → least)
PAIR_ORDER = ['both_expected_depth', 'one_expected_one_atypical',
              'both_atypical', 'asymmetric_coverage',
              'unresolved', 'concordant_half_depth_sum_valid']

PARALOG_TYPE_ORDER = ['l', 'dl', 'rl']
""")

add_code("""# Load v4 BUSCO-calibrated locus classification
loci_v4 = pd.read_csv(
    'read-depth-screen-v4-busco-calibrated/locus_classification_v4_busco.tsv', sep='\\t'
)
print(f'=== V4 BUSCO-calibrated Locus classification ===')
print(f'Total loci: {len(loci_v4)}  (parents: {(loci_v4["gene_type"] == "parent").sum()}, '
      f'paralogs: {(loci_v4["gene_type"] != "parent").sum()})')

# Separate paralogs and parents
paralog_loci = loci_v4[loci_v4['gene_type'] != 'parent'].copy()
parent_loci = loci_v4[loci_v4['gene_type'] == 'parent'].copy()

print(f'\\nLocus coverage categories (paralogs only, n={len(paralog_loci)}):')
for cat in LOCUS_ORDER:
    cnt = (paralog_loci['locus_coverage_category'] == cat).sum()
    if cnt > 0:
        pct = 100 * cnt / len(paralog_loci)
        print(f'  {cat:45s} {cnt:5d}  ({pct:5.1f}%)')

print(f'\\nThreshold version: {paralog_loci["threshold_version"].iloc[0]}')
print(f'BUSCO calibration N: {paralog_loci["busco_calibration_n"].iloc[0]}')
print(f'Expected depth: [{EXPECTED_LO}, {EXPECTED_HI}]')
print(f'Half-depth range: [{HALF_LO}, {HALF_HI}]')
print(f'Very-low range: [{VERY_LOW_LO}, {VERY_LOW_HI})')
""")

add_code("""# Load v4 BUSCO-calibrated pair integration
pairs_v4 = pd.read_csv(
    'read-depth-screen-v4-busco-calibrated/pair_integration_v4_busco.tsv', sep='\\t'
)
print(f'=== V4 BUSCO-calibrated Pair integration ===')
print(f'Total pairs: {len(pairs_v4)}')
for cat in PAIR_ORDER:
    cnt = (pairs_v4['pair_interpretation'] == cat).sum()
    if cnt > 0:
        pct = 100 * cnt / len(pairs_v4)
        print(f'  {cat:45s} {cnt:5d}  ({pct:5.1f}%)')
""")

add_code("""# Load genomic summary for reference (paralog_type already in locus table)
genomic_summary = pd.read_csv('paralog_genomic_summary.tsv', sep='\\t')
print(f'Genomic summary: {len(genomic_summary)} paralogs')
print(f'Paralog types in locus table: {paralog_loci[\"paralog_type\"].value_counts().to_dict()}')
""")

# ===========================================================================
# Panel A: Locus-level coverage classification criteria
# ===========================================================================
add_md("""## Panel A — Locus-level coverage classification criteria

Each of the 3,347 loci is classified independently by comparing its own
coverage (gene body + upstream/downstream flanks) to its chromosome-class
baseline (autosomal, chrX, or chrY).

**BUSCO-calibrated thresholds** (12,057 complete autosomal single-copy loci):

| Range | Interval | Rationale |
|-------|----------|-----------|
| Expected depth | [0.636, 1.394] | Central 95% of autosomal BUSCO |
| Half-depth-like | [0.30, 0.55] | ~C/2; false half-depth rate 0.6% among BUSCO |
| Very-low | [0.20, 0.30) | Below empirical BUSCO range (min = 0.394) |
| Intermediate | (0.55, 0.636) | Explicit gap between half-depth and expected |

**Decision logic** (from `calibrate_and_classify_v4.py`, first match wins):

1. Both flanks uninformative → `inconclusive`
2. Gene body unique < expected AND permissive > expected, both flanks normal → `unique_mapping_deficit`
3. Both flanks very-low AND combined very-low/half, both flanks informative → `very_low_coverage`
4. Both flanks half-depth AND combined half-depth, both flanks informative → `consistent_half_depth`
5. Both flanks expected AND combined expected, both flanks informative, gene body NOT half-depth → `normal_depth`
6. Any region at half-depth, very-low, or intermediate → `localized_or_asymmetric_low_depth`
7. Everything else → `inconclusive`
""")

add_code("""# ============================================================================
# Panel A — Locus coverage classification criteria (ax.table)
# ============================================================================

locus_criteria_examples = [
    (
        'EXPECTED DEPTH',
        '#2ECC40',
        'Both flanks informative and\\nboth within expected-depth\\nrange [0.636–1.394]',
        'Unique gene-body depth must\\nNOT fall in the half-depth\\nrange [0.30–0.55]',
        'Combined gene ±25 kb depth\\nmust be within expected-depth\\nrange [0.636–1.394]',
        LOCUS_INTERPRETATION['normal_depth'],
    ),
    (
        'HALF-DEPTH',
        '#E74C3C',
        'Both flanks informative and\\nboth within half-depth-like\\nrange [0.30–0.55]',
        'Gene body is not required to\\nbe half-depth; paralogous\\nsequence may distort gene-body\\nmapping',
        'Combined gene ±25 kb depth\\nmust also fall within half-depth\\nrange [0.30–0.55]',
        LOCUS_INTERPRETATION['consistent_half_depth'],
    ),
    (
        'VERY LOW COVERAGE',
        '#FF4136',
        'Both flanks informative and\\nboth within very-low range\\n[0.20–0.30)',
        'Gene body not required;\\ncombined region must also\\nbe very-low or half-depth',
        'Combined region also very-low\\nor half-depth-like relative\\nto chr-class baseline',
        LOCUS_INTERPRETATION['very_low_coverage'],
    ),
    (
        'LOCALIZED LOW DEPTH',
        '#FF851B',
        'At least one flank\\ninformative; one-sided,\\nasymmetric, intermediate,\\nor conflicting flank pattern',
        'Half-depth [0.30–0.55],\\nvery-low [0.20–0.30), or\\nintermediate (0.55–0.636)\\ngene-body depth may contribute',
        'Combined depth may be half-depth\\n[0.30–0.55], very-low\\n[0.20–0.30), intermediate\\n(0.55–0.636), or inconsistent\\nwith flanks',
        LOCUS_INTERPRETATION['localized_or_asymmetric_low_depth'],
    ),
    (
        'MAPPING AMBIGUITY',
        '#9B59B6',
        'Both flanks informative and\\nnear expected depth\\n[0.636–1.394]',
        'Unique gene-body depth < 0.636×C\\nwhile permissive gene-body depth\\n> 0.636×C — reads are present but\\ncannot be uniquely assigned',
        'Combined-region depth does not\\noverride the normal flank\\nsupport',
        LOCUS_INTERPRETATION['unique_mapping_deficit'],
    ),
    (
        'UNEVALUABLE',
        '#AAAAAA',
        'Both flanks uninformative\\n(≥50% zero-coverage) OR no\\ndefined pattern is matched',
        'Missing, intermediate, or\\ncontradictory evidence',
        'Missing, elevated (>1.394),\\nintermediate (0.55–0.636), or\\ncontradictory evidence',
        LOCUS_INTERPRETATION['inconclusive'],
    ),
]

col_labels = [
    'Category',
    'Flank requirements',
    'Gene-body role',
    'Combined-region role',
    'Interpretation',
]

# ============================================================================
# Build and draw locus criteria table with automatic wrapping
# ============================================================================

# Relative widths; must sum to 1.0
col_widths = [0.14, 0.18, 0.20, 0.18, 0.30]

# Approximate characters per wrapped line in each column.
wrap_widths = [24, 33, 38, 34, 52]


def normalize_text(text):
    return re.sub(r'\\s+', ' ', str(text)).strip()


def wrap_text(text, width):
    return textwrap.fill(
        normalize_text(text),
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
    )


# ---------------------------------------------------------------------------
# Print verification
# ---------------------------------------------------------------------------

print('=' * 100)
print('PANEL A — LOCUS CRITERIA TABLE')
print('=' * 100)

for i, (cat_name, cat_color, flank_req, gb_role, comb_role, interp) in enumerate(
        locus_criteria_examples, start=1):
    print(f'\\n{\"─\" * 100}')
    print(f'{i}. {normalize_text(cat_name)}')
    print(f'   Flank requirements:   {normalize_text(flank_req)}')
    print(f'   Gene-body role:       {normalize_text(gb_role)}')
    print(f'   Combined-region role: {normalize_text(comb_role)}')
    print(f'   Interpretation:       {normalize_text(interp)}')


# ---------------------------------------------------------------------------
# Build wrapped table contents
# ---------------------------------------------------------------------------

table_data = []
table_colors_locus = []
row_line_counts = []

for cat_name, cat_color, flank_req, gb_role, comb_role, interp in locus_criteria_examples:

    raw_values = [cat_name, flank_req, gb_role, comb_role, interp]

    wrapped_values = [
        wrap_text(value, width)
        for value, width in zip(raw_values, wrap_widths)
    ]

    table_data.append(wrapped_values)

    table_colors_locus.append([
        cat_color, '#FAFAFA', '#FAFAFA', '#FAFAFA', '#FAFAFA',
    ])

    max_lines = max(v.count('\\n') + 1 for v in wrapped_values)
    row_line_counts.append(max_lines)


# ---------------------------------------------------------------------------
# Create figure
# ---------------------------------------------------------------------------

fig, ax_table = plt.subplots(figsize=(20, 9))
ax_table.axis('off')

criteria_table = ax_table.table(
    cellText=table_data,
    colLabels=col_labels,
    cellColours=table_colors_locus,
    cellLoc='left',
    colLoc='center',
    loc='center',
    bbox=[0.015, 0.055, 0.97, 0.88],
)

criteria_table.auto_set_font_size(False)

n_body_rows = len(locus_criteria_examples)
n_total_rows = n_body_rows + 1
n_cols = len(col_labels)


# ---------------------------------------------------------------------------
# General cell formatting
# ---------------------------------------------------------------------------

for ci, width in enumerate(col_widths):
    for ri in range(n_total_rows):
        cell = criteria_table[ri, ci]
        cell.set_width(width)
        cell.set_edgecolor('#666666')
        cell.set_linewidth(0.6)
        cell.PAD = 0.065
        cell.get_text().set_wrap(True)
        cell.get_text().set_va('center')


# ---------------------------------------------------------------------------
# Header styling
# ---------------------------------------------------------------------------

for ci in range(n_cols):
    cell = criteria_table[0, ci]
    cell.set_facecolor('#333333')
    cell.set_text_props(
        color='white', fontweight='bold', fontsize=10,
        ha='center', va='center')
    cell.set_height(0.105)


# ---------------------------------------------------------------------------
# Body styling and dynamic row heights
# ---------------------------------------------------------------------------

height_per_line = 0.022
minimum_row_height = 0.105
maximum_row_height = 0.185

for ri in range(1, n_total_rows):
    cat_color = locus_criteria_examples[ri - 1][1]
    line_count = row_line_counts[ri - 1]

    row_height = min(
        maximum_row_height,
        max(minimum_row_height, 0.045 + height_per_line * line_count),
    )

    for ci in range(n_cols):
        criteria_table[ri, ci].set_height(row_height)

    # Category cell
    category_cell = criteria_table[ri, 0]
    category_cell.set_facecolor(cat_color)
    text_color = (
        '#222222' if cat_color in ('#FFDC00', '#F39C12') else 'white'
    )
    category_cell.set_text_props(
        color=text_color, fontweight='bold', fontsize=9,
        ha='left', va='center')

    # Remaining body cells
    for ci in range(1, n_cols):
        criteria_table[ri, ci].set_text_props(
            color='#222222', fontsize=8.4, ha='left', va='center')


# ---------------------------------------------------------------------------
# Title and export
# ---------------------------------------------------------------------------

ax_table.set_title(
    'BUSCO-calibrated HiFi read-depth criteria for paralog loci',
    fontsize=16, fontweight='bold', pad=18)

plt.savefig('panel_A_locus_criteria.png', dpi=300,
            bbox_inches='tight', facecolor='white')
plt.show()
print('\\nSaved: panel_A_locus_criteria.png')
""")

# ===========================================================================
# Panel B: Distribution of locus coverage categories (pie chart)
# ===========================================================================
add_md("""## Panel B — Distribution of locus coverage categories

Pie chart showing the proportion of annotated paralog loci assigned to
each coverage category. Each locus belongs to exactly one category.
""")

add_code("""# ============================================================================
# Panel B — Locus coverage category distribution (pie chart)
# ============================================================================

locus_counts = paralog_loci['locus_coverage_category'].value_counts()
locus_counts = locus_counts.reindex(
    [c for c in LOCUS_ORDER if c in locus_counts.index]
).fillna(0).astype(int)

fig, ax = plt.subplots(figsize=(7, 6))
colors_pie = [LOCUS_COLORS[c] for c in locus_counts.index]
wedges, texts, autotexts = ax.pie(
    locus_counts.values, labels=None, colors=colors_pie,
    autopct='%1.1f%%', startangle=90, pctdistance=0.55,
    textprops={'fontsize': 9})
for at in autotexts:
    at.set_fontweight('bold')

legend_labels = [f'{LOCUS_DISPLAY.get(c, c)} ({cnt})' for c, cnt in locus_counts.items()]
ax.legend(wedges, legend_labels, loc='center left',
          bbox_to_anchor=(1, 0.5), fontsize=7.5,
          title='Locus coverage category', title_fontsize=8.5)

plt.savefig('panel_B_locus_pie.png', dpi=300)
plt.show()
print('Saved: panel_B_locus_pie.png')

# Print counts
print('\\nLocus coverage categories (paralogs):')
print(f'{\"Category\":45s} {\"N\":>5s}  {\"%\":>6s}')
print('-' * 58)
for cat, cnt in locus_counts.items():
    print(f'{LOCUS_DISPLAY.get(cat, cat):45s} {cnt:5d}  {100*cnt/len(paralog_loci):5.1f}%')
""")

# ===========================================================================
# Panel C: Coverage category by paralog type (stacked bar %)
# ===========================================================================
add_md("""## Panel C — Coverage category by paralog type

Stacked bar chart showing the distribution of locus coverage categories
across paralog classes (l: Liftoff-only, dl: de novo + Liftoff,
rl: reference-based Liftoff).
""")

add_code("""# ============================================================================
# Panel C — Locus coverage category by paralog type (stacked bar %)
# ============================================================================

paralog_only = paralog_loci[paralog_loci['paralog_type'].notna() & (paralog_loci['paralog_type'] != '')].copy()

ct_abs = pd.crosstab(
    paralog_only['locus_coverage_category'],
    paralog_only['paralog_type']
)
ct_abs = ct_abs.reindex(
    index=[c for c in LOCUS_ORDER if c in ct_abs.index],
    columns=PARALOG_TYPE_ORDER
).fillna(0).astype(int)
ct_pct = ct_abs.div(ct_abs.sum(axis=0), axis=1) * 100

fig, ax = plt.subplots(figsize=(7, 5.5))
bottom = np.zeros(len(PARALOG_TYPE_ORDER))
for cat in ct_pct.index:
    vals = ct_pct.loc[cat].values
    ax.bar(range(len(PARALOG_TYPE_ORDER)), vals, bottom=bottom,
           color=LOCUS_COLORS[cat], label=LOCUS_DISPLAY.get(cat, cat),
           width=0.6, edgecolor='white', linewidth=0.5)
    for i, v in enumerate(vals):
        if v > 5:
            ax.text(i, bottom[i] + v/2, f'{v:.0f}%', ha='center',
                   va='center', fontsize=8, fontweight='bold')
    bottom += vals

ax.set_xticks(range(len(PARALOG_TYPE_ORDER)))
ax.set_xticklabels([f'{t}\\n(n={int(ct_abs[t].sum())})' for t in PARALOG_TYPE_ORDER])
ax.set_ylabel('Percentage of paralog loci')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=7,
          title='Locus category', title_fontsize=8)
ax.set_ylim(0, 105)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('panel_C_locus_by_paralog_type.png', dpi=300)
plt.show()
print('Saved: panel_C_locus_by_paralog_type.png')

# Print cross-tab
ct_print = ct_abs.copy()
ct_print['Total'] = ct_print.sum(axis=1)
ct_print.loc['Total'] = ct_print.sum(axis=0)
print('\\n=== Locus coverage category x Paralog type ===')
print(ct_print.to_string())
""")

# ===========================================================================
# Panels D-I: Representative locus examples
# ===========================================================================
add_md("""## Panels D–I — Representative locus examples

One representative paralog locus for each coverage category. Each panel
shows ONLY the paralog locus (not the parent). Coverage values are shown
for four regions: upstream flank, gene body, downstream flank, and combined
gene ±25 kb. BUSCO-derived thresholds are overlaid as horizontal bands.
""")

add_code("""# ============================================================================
# Panels D-I — Representative locus examples (one per category)
# ============================================================================

# Select one representative paralog per category
np.random.seed(42)
locus_examples = {}
for cat in LOCUS_ORDER:
    subset = paralog_loci[paralog_loci['locus_coverage_category'] == cat]
    if len(subset) == 0:
        continue
    # Prefer 'dl' type for examples (most common paralog type)
    dl_subset = subset[subset['paralog_type'] == 'dl']
    if len(dl_subset) > 0:
        pick = dl_subset.sample(1, random_state=42).iloc[0]
    else:
        pick = subset.sample(1, random_state=42).iloc[0]
    locus_examples[cat] = pick
    print(f'{cat:45s} → {pick[\"gene_name\"]:25s}  type={pick.get(\"paralog_type\", \"?\")}')

# Panel letters for each category
panel_letters = {
    'normal_depth': 'D',
    'consistent_half_depth': 'E',
    'very_low_coverage': 'F',
    'localized_or_asymmetric_low_depth': 'G',
    'unique_mapping_deficit': 'H',
    'inconclusive': 'I',
}

# Region labels for x-axis
region_labels = ['Upstream\\nflank', 'Gene\\nbody', 'Downstream\\nflank', 'Combined\\n±25 kb']

for cat in LOCUS_ORDER:
    if cat not in locus_examples:
        continue

    row = locus_examples[cat]
    gene_name = row['gene_name']
    panel = panel_letters.get(cat, '?')

    # Coverage values (unique and permissive)
    up_uq   = row.get('upstream_unique_chr_norm', np.nan)
    gb_uq   = row.get('gene_body_unique_chr_norm', np.nan)
    dn_uq   = row.get('downstream_unique_chr_norm', np.nan)
    comb_uq = row.get('combined_unique_chr_norm', np.nan)
    up_perm   = row.get('upstream_permissive_chr_norm', np.nan)
    gb_perm   = row.get('gene_body_permissive_chr_norm', np.nan)
    dn_perm   = row.get('downstream_permissive_chr_norm', np.nan)
    comb_perm = row.get('combined_permissive_chr_norm', np.nan)

    unique_vals = [up_uq, gb_uq, dn_uq, comb_uq]
    permissive_vals = [up_perm, gb_perm, dn_perm, comb_perm]

    # Diagnostic flags
    up_half = row.get('upstream_is_half', False)
    dn_half = row.get('downstream_is_half', False)
    up_dipl = row.get('upstream_is_diploid', False)
    dn_dipl = row.get('downstream_is_diploid', False)
    up_vlow = row.get('upstream_is_very_low', False)
    dn_vlow = row.get('downstream_is_very_low', False)
    sym_half = row.get('symmetric_half_depth', False)
    map_def = row.get('has_mapping_deficit', False)
    both_info = row.get('both_flanks_informative', False)
    comb_half = row.get('combined_is_half', False)
    comb_dipl = row.get('combined_is_diploid', False)
    gb_dipl = row.get('gene_body_is_diploid', False)

    fig, ax = plt.subplots(figsize=(6, 4.5))

    x = np.arange(len(region_labels))
    width = 0.30

    bars_uq = ax.bar(x - width/2, unique_vals, width,
                     color='#3498DB', label='Unique (MAPQ≥20)',
                     edgecolor='white', linewidth=0.5)
    bars_perm = ax.bar(x + width/2, permissive_vals, width,
                       color='#95A5A6', label='Permissive (MAPQ≥0)',
                       edgecolor='white', linewidth=0.5)

    # Value labels on bars
    for i, (vu, vp) in enumerate(zip(unique_vals, permissive_vals)):
        if np.isfinite(vu):
            ax.text(i - width/2, vu + 0.03, f'{vu:.2f}', ha='center',
                   va='bottom', fontsize=7, fontweight='bold', color='#1a5276')
        if np.isfinite(vp):
            ax.text(i + width/2, vp + 0.03, f'{vp:.2f}', ha='center',
                   va='bottom', fontsize=7, color='#5d6d7e')

    # Threshold bands
    ax.axhspan(EXPECTED_LO, EXPECTED_HI, alpha=0.10, color='#2ECC40', zorder=0)
    ax.axhspan(HALF_LO, HALF_HI, alpha=0.10, color='#E74C3C', zorder=0)
    ax.axhspan(VERY_LOW_LO, VERY_LOW_HI, alpha=0.10, color='#FF4136', zorder=0)

    # Threshold lines
    for thresh, ls, color, label in [
        (EXPECTED_LO, '--', '#27AE60', f'Expected lo ({EXPECTED_LO})'),
        (EXPECTED_HI, '--', '#27AE60', f'Expected hi ({EXPECTED_HI})'),
        (HALF_LO, ':', '#C0392B', f'Half-depth lo ({HALF_LO})'),
        (HALF_HI, ':', '#C0392B', f'Half-depth hi ({HALF_HI})'),
        (VERY_LOW_LO, '-.', '#E74C3C', f'Very-low ({VERY_LOW_LO})'),
    ]:
        ax.axhline(y=thresh, color=color, linestyle=ls, linewidth=0.8, alpha=0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(region_labels, fontsize=8)
    ax.set_ylabel('Normalized coverage (× chr-class baseline)', fontsize=9)
    ax.set_ylim(0, max(max(unique_vals), max(permissive_vals), EXPECTED_HI) * 1.25 + 0.1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend with threshold explanation
    legend_elements = [
        mpatches.Patch(facecolor='#3498DB', label='Unique (MAPQ≥20)'),
        mpatches.Patch(facecolor='#95A5A6', label='Permissive (MAPQ≥0)'),
        mpatches.Patch(facecolor='#2ECC40', alpha=0.2, label=f'Expected [{EXPECTED_LO}–{EXPECTED_HI}]'),
        mpatches.Patch(facecolor='#E74C3C', alpha=0.2, label=f'Half-depth [{HALF_LO}–{HALF_HI}]'),
        mpatches.Patch(facecolor='#FF4136', alpha=0.2, label=f'Very-low [{VERY_LOW_LO}–{VERY_LOW_HI})'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1),
             fontsize=6.5, title=f'{gene_name}', title_fontsize=8)

    # Classification summary as text box
    flag_text = []
    if both_info:
        flag_text.append('flanks informative')
    else:
        flag_text.append('flanks NOT informative')
    if sym_half:
        flag_text.append('symmetric half-depth')
    if map_def:
        flag_text.append('mapping deficit')
    if comb_half:
        flag_text.append('combined half')
    if comb_dipl:
        flag_text.append('combined expected')
    if gb_dipl:
        flag_text.append('gene body expected')
    if up_half or dn_half:
        sides = []
        if up_half: sides.append('up')
        if dn_half: sides.append('down')
        flag_text.append(f'flank half: {\",\".join(sides)}')
    if up_dipl or dn_dipl:
        sides = []
        if up_dipl: sides.append('up')
        if dn_dipl: sides.append('down')
        flag_text.append(f'flank expected: {\",\".join(sides)}')
    if up_vlow or dn_vlow:
        sides = []
        if up_vlow: sides.append('up')
        if dn_vlow: sides.append('down')
        flag_text.append(f'flank very-low: {\",\".join(sides)}')

    cat_display = LOCUS_DISPLAY.get(cat, cat)
    ax.text(0.02, 0.98, f'{cat_display}\\n' + '\\n'.join(flag_text),
            transform=ax.transAxes, fontsize=6.5, verticalalignment='top',
            fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    fname = f'panel_{panel}_{cat}.png'
    plt.savefig(fname, dpi=300)
    plt.show()
    print(f'Saved: {fname}')
""")

# ===========================================================================
# Panel J: Pair-level interpretation framework
# ===========================================================================
add_md("""## Panel J — Pair-level interpretation framework

After each locus is independently classified, parent–paralog pairs are
assigned an interpretation by combining the two locus classifications.
The sum check (parent + paralog combined depth ≈ 1.0) is used as a
secondary validator: two half-depth loci that sum to the expected
single-copy baseline are consistent with a single diploid locus
represented as two haplotypic sequences.

This panel explains how two independently classified loci are combined
into biological interpretations. It is intentionally an interpretation
framework rather than another coverage classifier.
""")

add_code("""# ============================================================================
# Panel J — Pair-level interpretation framework (ax.table)
# ============================================================================

# Ordered by count (most → least frequent)
# Tuple: (display_name, color, locus_combination, biological_interpretation)
# Sum-check requirements are rolled into the locus-combination column.
pair_criteria_examples = [
    (
        'BOTH LOCI AT EXPECTED DEPTH',
        '#2ECC40',
        'Both loci = EXPECTED DEPTH',
        PAIR_INTERPRETATION_TEXT['both_expected_depth'],
    ),
    (
        'ONE EXPECTED, ONE ATYPICAL',
        '#FFDC00',
        'One locus = EXPECTED DEPTH,\\nother = LOCALIZED LOW DEPTH,\\nUNEVALUABLE, or MAPPING AMBIGUITY',
        PAIR_INTERPRETATION_TEXT['one_expected_one_atypical'],
    ),
    (
        'BOTH ATYPICAL',
        '#F39C12',
        'At least one locus = LOCALIZED LOW DEPTH,\\nVERY LOW COVERAGE, or MAPPING AMBIGUITY;\\nthe other may be UNEVALUABLE;\\nnot both EXPECTED DEPTH or both HALF-DEPTH',
        PAIR_INTERPRETATION_TEXT['both_atypical'],
    ),
    (
        'ASYMMETRIC COVERAGE',
        '#FF851B',
        'One locus = EXPECTED DEPTH,\\nother = HALF-DEPTH\\nor VERY LOW COVERAGE',
        PAIR_INTERPRETATION_TEXT['asymmetric_coverage'],
    ),
    (
        'INDETERMINATE',
        '#AAAAAA',
        'Both loci = UNEVALUABLE',
        PAIR_INTERPRETATION_TEXT['unresolved'],
    ),
    (
        'CONCORDANT HALF-DEPTH',
        '#E74C3C',
        'Both loci = HALF-DEPTH\\nSum check ≈ 1.0 (same chr class,\\nnon-overlapping, both flanks\\ninformative)',
        PAIR_INTERPRETATION_TEXT['concordant_half_depth_sum_valid'],
    ),
]

col_labels_pair = [
    'Interpretation',
    'Locus combination',
    'Biological interpretation',
]

# ============================================================================
# Build and draw pair-level criteria table with automatic text wrapping
# ============================================================================

# Relative column widths; must sum to 1.0
col_widths_pair = [0.19, 0.33, 0.48]

# Approximate characters per line for each column
wrap_widths_pair = [30, 52, 72]


def normalize_text(text):
    return re.sub(r'\\s+', ' ', str(text)).strip()


def wrap_text(text, width):
    return textwrap.fill(
        normalize_text(text),
        width=width,
        break_long_words=False,
        break_on_hyphens=False,
    )


# ---------------------------------------------------------------------------
# Print verification
# ---------------------------------------------------------------------------

print('=' * 100)
print('PANEL J — PAIR INTERPRETATION FRAMEWORK')
print('=' * 100)

for i, (interp_name, interp_color, locus_combo, bio_interp) in enumerate(
        pair_criteria_examples, start=1):
    print(f'\\n{\"─\" * 100}')
    print(f'{i}. {normalize_text(interp_name)}')
    print(f'   Locus combination:       {normalize_text(locus_combo)}')
    print(f'   Biological interpretation: {normalize_text(bio_interp)}')


# ---------------------------------------------------------------------------
# Build automatically wrapped table contents
# ---------------------------------------------------------------------------

table_data_pair = []
table_colors_pair = []
row_line_counts = []

for interp_name, interp_color, locus_combo, bio_interp in pair_criteria_examples:

    raw_values = [interp_name, locus_combo, bio_interp]

    wrapped_values = [
        wrap_text(value, width)
        for value, width in zip(raw_values, wrap_widths_pair)
    ]

    table_data_pair.append(wrapped_values)
    table_colors_pair.append([interp_color, '#FAFAFA', '#FAFAFA'])

    max_lines = max(v.count('\\n') + 1 for v in wrapped_values)
    row_line_counts.append(max_lines)


# ---------------------------------------------------------------------------
# Create figure
# ---------------------------------------------------------------------------

fig, ax_table = plt.subplots(figsize=(18, 7.5))
ax_table.axis('off')

criteria_table = ax_table.table(
    cellText=table_data_pair,
    colLabels=col_labels_pair,
    cellColours=table_colors_pair,
    cellLoc='left',
    colLoc='center',
    loc='center',
    bbox=[0.02, 0.06, 0.96, 0.88],
)

criteria_table.auto_set_font_size(False)

n_body_rows = len(pair_criteria_examples)
n_total_rows = n_body_rows + 1
n_cols = len(col_labels_pair)


# ---------------------------------------------------------------------------
# General cell formatting
# ---------------------------------------------------------------------------

for ci, width in enumerate(col_widths_pair):
    for ri in range(n_total_rows):
        cell = criteria_table[ri, ci]
        cell.set_width(width)
        cell.set_edgecolor('#666666')
        cell.set_linewidth(0.6)
        cell.PAD = 0.07
        cell.get_text().set_wrap(True)
        cell.get_text().set_va('center')


# ---------------------------------------------------------------------------
# Header styling
# ---------------------------------------------------------------------------

for ci in range(n_cols):
    cell = criteria_table[0, ci]
    cell.set_facecolor('#333333')
    cell.set_text_props(
        color='white', fontweight='bold', fontsize=10,
        ha='center', va='center')
    cell.set_height(0.105)


# ---------------------------------------------------------------------------
# Body styling and dynamic row heights
# ---------------------------------------------------------------------------

height_per_line = 0.026
minimum_row_height = 0.105
maximum_row_height = 0.175

for ri in range(1, n_total_rows):
    interp_color = pair_criteria_examples[ri - 1][1]
    line_count = row_line_counts[ri - 1]

    row_height = min(
        maximum_row_height,
        max(minimum_row_height, 0.050 + height_per_line * line_count),
    )

    for ci in range(n_cols):
        criteria_table[ri, ci].set_height(row_height)

    # Colored interpretation cell
    category_cell = criteria_table[ri, 0]
    category_cell.set_facecolor(interp_color)
    text_color = (
        '#222222' if interp_color in ('#FFDC00', '#F39C12') else 'white'
    )
    category_cell.set_text_props(
        color=text_color, fontweight='bold', fontsize=9,
        ha='left', va='center')

    # Other body cells
    for ci in range(1, n_cols):
        criteria_table[ri, ci].set_text_props(
            color='#222222', fontsize=8.7, ha='left', va='center')


# ---------------------------------------------------------------------------
# Title and export
# ---------------------------------------------------------------------------

ax_table.set_title(
    'Pair-level integration of locus coverage classifications',
    fontsize=15, fontweight='bold', pad=16)

plt.savefig('panel_J_pair_criteria.png', dpi=300,
            bbox_inches='tight', facecolor='white')
plt.show()
print('\\nSaved: panel_J_pair_criteria.png')
""")

# ===========================================================================
# Assemble notebook
# ===========================================================================
add_md("""## Output files

| Panel | File | Description |
|-------|------|-------------|
| A | `panel_A_locus_criteria.png` | Locus-level coverage classification criteria table |
| B | `panel_B_locus_pie.png` | Distribution of locus coverage categories (pie chart) |
| C | `panel_C_locus_by_paralog_type.png` | Coverage category by paralog type (stacked bar %) |
| D | `panel_D_normal_depth.png` | Representative example: EXPECTED DEPTH |
| E | `panel_E_consistent_half_depth.png` | Representative example: HALF-DEPTH |
| F | `panel_F_very_low_coverage.png` | Representative example: VERY LOW COVERAGE |
| G | `panel_G_localized_or_asymmetric_low_depth.png` | Representative example: LOCALIZED LOW DEPTH |
| H | `panel_H_unique_mapping_deficit.png` | Representative example: MAPPING AMBIGUITY |
| I | `panel_I_inconclusive.png` | Representative example: UNEVALUABLE |
| J | `panel_J_pair_criteria.png` | Pair-level interpretation framework table |
""")

add_code("""print("\\n" + "="*60)
print("ALL PANELS GENERATED")
print("="*60)
print("\\nPanel files:")
import glob
for f in sorted(glob.glob('panel_*.png')):
    print(f'  {f}')
""")

# ===========================================================================
# Assemble
# ===========================================================================
nb.cells = cells

nbf.write(nb, 'category_exploration_v4_busco.ipynb')
print(f'Notebook written with {len(cells)} cells')
