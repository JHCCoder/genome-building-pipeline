#!/usr/bin/env python3
"""
Heatmap of all paralog expression (TPM) across 29 RNA-seq samples.
- Rows: all 2,252 paralogs, clustered within family blocks
- Columns: all 29 samples, grouped by tissue
- Row labels hidden (too many), family boundaries marked
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
counts_raw = pd.read_csv('paralog_families_counts.tsv', sep='\t', comment='#', index_col=0)
sample_cols_orig = [c for c in counts_raw.columns if c.startswith('SRR')]
counts = counts_raw[sample_cols_orig].copy()
counts.columns = [c.replace('_Aligned.sortedByCoord.out.bam', '') for c in counts.columns]
sample_list = list(counts.columns)

families_df = pd.read_csv('paralog_families.tsv', sep='\t').set_index('gene_id')

# Merge and filter to paralogs only
meta_cols = ['family', 'gene_name', 'gene_type', 'paralog_type', 'length']
merged = counts.join(families_df[meta_cols], how='inner')
paralogs = merged[merged['gene_type'] != 'parent'].copy()
print(f"Paralogs: {len(paralogs)}")

# TPM normalization
gene_lengths = paralogs['length']
rpk = paralogs[sample_list].div(gene_lengths / 1000, axis=0)
tpm = rpk.div(rpk.sum(axis=0), axis=1) * 1e6

# Tissue metadata
sample_tissue = {
    'SRR17216293': 'skin',    'SRR17216294': 'brain',   'SRR17216295': 'heart',
    'SRR17216296': 'kidney',  'SRR17216297': 'liver',   'SRR17216298': 'lung',
    'SRR17216299': 'skin',    'SRR17216300': 'brain',   'SRR17216301': 'heart',
    'SRR17216302': 'kidney',  'SRR17216303': 'liver',   'SRR17216304': 'lung',
    'SRR17216305': 'skin',    'SRR17216306': 'brain',   'SRR17216307': 'heart',
    'SRR17216308': 'kidney',  'SRR17216309': 'liver',   'SRR17216310': 'lung',
    'SRR17216311': 'brain',   'SRR17216312': 'heart',   'SRR17216313': 'kidney',
    'SRR17216314': 'kidney',  'SRR17216315': 'lung',    'SRR17216316': 'skin',
    'SRR17216317': 'heart',   'SRR17216318': 'kidney',  'SRR17216319': 'liver',
    'SRR17216320': 'brain',   'SRR17216321': 'lung'
}
tissues = ['skin', 'brain', 'heart', 'kidney', 'liver', 'lung']
tissue_colors = {
    'skin': '#E8A87C', 'brain': '#C38D9E', 'heart': '#E27D60',
    'kidney': '#85CDCA', 'liver': '#D4A574', 'lung': '#41B3A3'
}

# Sort samples by tissue
sample_order = sorted(sample_list, key=lambda s: (tissues.index(sample_tissue[s]), s))
tpm_sorted = tpm[sample_order]

# ---------------------------------------------------------------------------
# Row ordering: cluster families, keep family members together
# ---------------------------------------------------------------------------
# Compute per-family mean expression vector across samples
family_mean = tpm.groupby(paralogs['family']).mean()
# Cluster families
family_linkage = linkage(family_mean, method='ward', metric='euclidean')
family_dendro = dendrogram(family_linkage, no_plot=True)
family_order = [family_mean.index[i] for i in family_dendro['leaves']]

# Within each family, sort paralogs by total expression
row_order = []
family_boundaries = []  # (start_idx, end_idx, family_name)
for fam in family_order:
    fam_mask = paralogs['family'] == fam
    fam_paralogs = tpm_sorted[fam_mask]
    # Sort paralogs within family by mean expression
    fam_order = fam_paralogs.mean(axis=1).sort_values(ascending=False).index
    start = len(row_order)
    row_order.extend(fam_order)
    end = len(row_order)
    family_boundaries.append((start, end, fam))

tpm_ordered = tpm_sorted.loc[row_order]

# Log transform for visualization
log_tpm = np.log2(tpm_ordered.values + 0.1)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(16, 30))

# Main heatmap
ax_heatmap = fig.add_axes([0.08, 0.05, 0.82, 0.93])

im = ax_heatmap.imshow(log_tpm, aspect='auto', cmap='YlOrRd',
                        vmin=-3, vmax=np.percentile(log_tpm, 99),
                        interpolation='none')

# Column labels: sample names colored by tissue
ax_heatmap.set_xticks(range(len(sample_order)))
ax_heatmap.set_xticklabels(sample_order, rotation=90, fontsize=5, ha='center')
for i, s in enumerate(sample_order):
    ax_heatmap.get_xticklabels()[i].set_color(tissue_colors[sample_tissue[s]])

# Tissue group labels on top
prev_tissue = None
boundaries = []
for i, s in enumerate(sample_order):
    t = sample_tissue[s]
    if t != prev_tissue:
        boundaries.append(i)
        prev_tissue = t
boundaries.append(len(sample_order))

for i in range(len(boundaries) - 1):
    mid = (boundaries[i] + boundaries[i+1]) / 2
    tissue_name = sample_tissue[sample_order[boundaries[i]]]
    ax_heatmap.text(mid, -1.8, tissue_name, ha='center', va='bottom',
                    fontsize=8, fontweight='bold', color=tissue_colors[tissue_name],
                    transform=ax_heatmap.transData)

# Draw tissue separators
for b in boundaries[1:-1]:
    ax_heatmap.axvline(x=b - 0.5, color='white', linewidth=1.5, alpha=0.8)

# Y-axis: hide individual row labels, show family labels
ax_heatmap.set_yticks([])
# Add family labels on the right side
for start, end, fam_name in family_boundaries:
    mid = (start + end) / 2
    # Only label families with >= 2 members or every Nth family to avoid crowding
    if end - start >= 3:
        ax_heatmap.text(len(sample_order) + 0.5, mid, fam_name,
                        fontsize=3.5, va='center', ha='left', color='#333333')

# Color bar
cbar_ax = fig.add_axes([0.91, 0.05, 0.015, 0.93])
cbar = plt.colorbar(im, cax=cbar_ax)
cbar.set_label('log2(TPM + 0.1)', fontsize=9)

# Title
fig.suptitle(f'Paralog Expression Across Tissues\n'
             f'{len(paralogs)} paralogs, {len(family_order)} families, 29 samples',
             fontsize=12, fontweight='bold', y=0.98)

# Tissue legend at bottom
legend_handles = [mpatches.Patch(color=tissue_colors[t], label=t) for t in tissues]
fig.legend(handles=legend_handles, loc='lower center', ncol=6, fontsize=8,
           frameon=False, bbox_to_anchor=(0.5, 0.01))

plt.savefig('paralog_expression_heatmap_all.png', dpi=300)
plt.savefig('paralog_expression_heatmap_all.pdf', dpi=300)
print("Saved: paralog_expression_heatmap_all.png")
print("Saved: paralog_expression_heatmap_all.pdf")
