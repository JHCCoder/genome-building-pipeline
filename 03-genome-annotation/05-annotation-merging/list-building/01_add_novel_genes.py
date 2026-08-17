# %% [markdown]
# # Scenario 1 — Add novel de novo genes
#
# Adds functionally-annotated **de novo (Braker3) genes** that are absent from the
# Liftoff-transferred annotation. A de novo gene is added when it satisfies all
# three criteria:
#
# 1. **No CDS overlap** with any Liftoff gene (strand-aware, same-strand overlap).
# 2. **CDS length ≥ 60 %** of the same-named **mouse (mm39 / GRCm39)** homolog.
# 3. **Expression > 0** — mean RNA-seq signal over the gene's CDS (BigWig).
#
# A final de-duplication step drops a de novo gene if it overlaps a *same-named,
# same-strand* Liftoff gene at the gene level. Paralogs are then suffixed
# `-dl1`, `-dl2`, … so every gene name is unique.
#
# **Output:** `gene_gBraker_toAdd_unique.tsv` (`gene_id<TAB>new_gene_name`).
#
# **Run after** `02_rename_loc_genes` — its rename list is used to name the
# Liftoff genes before the same-name ambiguity filter (see step 2 below).

# %% [markdown]
# ## 0. Configuration
#
# Edit the paths below for your run. The defaults reproduce the *Octodon degus*
# annotation.

# %%
# --- inputs ---
LIFTOFF_GFF_GENES = "hifiasm-041425-scaffolded-chrAssigned-mito.gff"  # Liftoff gene table (naming)
LIFTOFF_GFF_CDS   = "gffl1.gff"                                        # Liftoff GFF used for CDS overlap
DENOVO_GFF        = "gffd2_sorted.gff"                                 # Braker3 de novo (functionally named)
RENAME_LIST       = "gene_LOC_nameChange_unique.tsv"                   # produced by 02_rename_loc_genes
MOUSE_GFF         = "GCF_000001635.27_GRCm39_genomic.gff"              # mm39 / GRCm39 (ortholog CDS lengths)
BIGWIG            = "merged.bw"                                        # RNA-seq expression track
MIN_MOUSE_CDS_RATIO = 0.6                                              # criterion 2 threshold

# --- outputs ---
OUT_ADD_LIST = "gene_gBraker_toAdd_unique.tsv"

# %% [markdown]
# ## 1. Imports and helper functions

# %%
import pandas as pd
import numpy as np
import pyBigWig
from collections import defaultdict
from intervaltree import IntervalTree


def parse_sorted_gff(gff_file):
    """Fast parser for a *sorted* GFF3 file.

    Assumes the file is sorted so that a gene is followed by its mRNAs, and each
    mRNA is followed by its CDS/exon features.
    """
    gene_dict = {}
    current_gene = None
    current_mrna = None

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields

            attr_dict = {}
            for attr in filter(None, attributes.split(";")):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attr_dict[key] = val
            if "ID" not in attr_dict:
                continue

            if feature_type == "gene":
                current_gene = attr_dict["ID"]
                gene_dict[current_gene] = {
                    "type": "gene", "seqid": seqid, "source": source,
                    "start": int(start), "end": int(end), "score": score,
                    "strand": strand, "phase": phase, "attributes": attr_dict,
                    "children": {},
                }
                current_mrna = None
            elif feature_type == "mRNA":
                if current_gene is None:
                    continue
                current_mrna = attr_dict["ID"]
                gene_dict[current_gene]["children"][current_mrna] = {
                    "type": "mRNA", "seqid": seqid, "source": source,
                    "start": int(start), "end": int(end), "score": score,
                    "strand": strand, "phase": phase, "attributes": attr_dict,
                    "children": {},
                }
            elif feature_type in ("CDS", "exon", "start_codon", "stop_codon"):
                if current_mrna is None:
                    continue
                feature_id = attr_dict["ID"]
                feature_data = {
                    "type": feature_type, "seqid": seqid, "source": source,
                    "start": int(start), "end": int(end), "score": score,
                    "strand": strand, "phase": phase, "attributes": attr_dict,
                }
                # CDS IDs are made unique by appending their coordinates.
                if feature_type == "CDS":
                    feature_id = f"{attr_dict['ID']}|{start}-{end}"
                gene_dict[current_gene]["children"][current_mrna]["children"][feature_id] = feature_data

    return gene_dict


def extract_cds(gff_dict, source_gff):
    """Flatten CDS features of a parsed GFF into a DataFrame."""
    cds_features = []
    for gene_id, gene_data in gff_dict.items():
        for mrna_id, mrna_data in gene_data["children"].items():
            mrna_attrs = mrna_data.get("attributes", {})
            protein_name = mrna_attrs.get("Name", "")
            product_name = mrna_attrs.get("product", "")
            for feature_id, feature_data in mrna_data["children"].items():
                if feature_data["type"] == "CDS":
                    cds_features.append({
                        "gene_id": gene_id,
                        "gene_name": gene_data["attributes"].get("Name", ""),
                        "mrna_id": mrna_id,
                        "mrna_name": protein_name,
                        "protein_product": product_name,
                        "cds_id": feature_id,
                        "chr": feature_data["seqid"],
                        "start": int(feature_data["start"]),
                        "end": int(feature_data["end"]),
                        "strand": feature_data["strand"],
                        "phase": feature_data.get("phase", "."),
                        "source_gff": source_gff,
                    })
    return pd.DataFrame(cds_features)


def find_non_overlapping_genes_with_cds(gff1_dict, gff2_dict):
    """Return CDS of de novo genes (gff1) whose CDS never overlap any Liftoff
    (gff2) CDS on the same chromosome and strand."""
    cds1 = extract_cds(gff1_dict, "gff1")
    cds2 = extract_cds(gff2_dict, "gff2")

    chr_strand_trees = {}
    for _, row in cds2.iterrows():
        key = (row["chr"], row["strand"])
        chr_strand_trees.setdefault(key, IntervalTree()).addi(row["start"], row["end"] + 1)

    genes_with_overlaps = set()
    for _, row in cds1.iterrows():
        key = (row["chr"], row["strand"])
        if key in chr_strand_trees and chr_strand_trees[key].overlaps(row["start"], row["end"] + 1):
            genes_with_overlaps.add(row["gene_id"])

    non_overlapping = cds1[~cds1["gene_id"].isin(genes_with_overlaps)]
    return non_overlapping.reset_index(drop=True)


def filter_entire_genes_with_overlaps(cds_df, gff_df):
    """Drop de novo CDS whose gene overlaps a *same-named, same-strand* Liftoff
    gene at the gene level (removes duplicates of genes already present)."""
    def create_keys(df, chr_col="chr", strand_col="strand", name_col="gene_name"):
        return list(zip(df[chr_col], df[strand_col], df[name_col]))

    cds_keys = create_keys(cds_df)
    gff_keys = create_keys(gff_df, chr_col="seqid", strand_col="strand", name_col="gene_name")

    gene_trees = {}
    for key, (_, gene_row) in zip(gff_keys, gff_df.iterrows()):
        gene_trees.setdefault(key, IntervalTree()).addi(gene_row["start"], gene_row["end"] + 1)

    genes_to_exclude = set()
    for key, (_, cds_row) in zip(cds_keys, cds_df.iterrows()):
        if key in gene_trees and gene_trees[key].overlaps(cds_row["start"], cds_row["end"] + 1):
            genes_to_exclude.add(key)

    mask = [k not in genes_to_exclude for k in cds_keys]
    return cds_df[mask].copy().reset_index(drop=True)


def get_bw_stats(bw, chrom, start, end, stat="mean"):
    """Mean BigWig signal over an interval (or NaN if the region is unmapped)."""
    try:
        values = bw.values(chrom, start, end)
        values = [v for v in values if v is not None]
        if not values:
            return 0.0
        return float(np.mean(values))
    except Exception:
        return np.nan


def calculate_gene_cds_lengths(df):
    """Total CDS length per gene (sum of CDS row spans)."""
    return df.groupby("gene_id").apply(lambda x: (x["end"] - x["start"] + 1).sum()).to_dict()


def parse_mouse_gff(mouse_gff_path):
    """Total CDS length per mouse gene (uppercase keys) from a mouse GFF3."""
    mouse_genes = defaultdict(list)
    with open(mouse_gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = {}
            for item in parts[8].replace("%3B", ";").split(";"):
                if "=" in item:
                    key, val = [x.strip() for x in item.split("=", 1)]
                    attrs[key.upper()] = val.upper()
            gene_name = next((attrs[f] for f in ("GENE", "NAME", "GENE_NAME", "GENE_SYMBOL") if f in attrs), None)
            if gene_name:
                mouse_genes[gene_name].append(int(parts[4]) - int(parts[3]) + 1)
    return {gene: sum(lengths) for gene, lengths in mouse_genes.items()}


def filter_genes_with_exons(your_genes_df, mouse_cds_dict, threshold=0.6):
    """Score each de novo gene against its same-named mouse homolog.

    ``status`` is ``"complete"`` when the de novo CDS length is ≥ ``threshold``
    of the mouse homolog's CDS length, ``"incomplete"`` otherwise, and
    ``"no_ortholog"`` when no same-named mouse gene is found.
    """
    your_gene_cds = calculate_gene_cds_lengths(your_genes_df)
    gene_meta = your_genes_df.drop_duplicates("gene_id").set_index("gene_id")
    gene_meta["GENE_NAME_UPPER"] = gene_meta["gene_name"].str.upper()
    gene_meta["MRNA_NAME_UPPER"] = gene_meta["mrna_name"].str.upper()

    results = []
    for gene_id, total_cds in your_gene_cds.items():
        row = gene_meta.loc[gene_id]
        match_keys = [row["GENE_NAME_UPPER"], row["MRNA_NAME_UPPER"],
                      row["gene_name"].upper(), row["mrna_name"].upper()]
        mouse_total = next((mouse_cds_dict[k] for k in match_keys if k in mouse_cds_dict), None)
        if mouse_total:
            ratio = total_cds / mouse_total
            status = "complete" if ratio >= threshold else "incomplete"
        else:
            ratio = np.nan
            status = "no_ortholog"
        results.append({
            "gene_id": gene_id, "gene_name": row["gene_name"], "mrna_name": row["mrna_name"],
            "your_total_cds": total_cds, "mouse_total_cds": mouse_total,
            "length_ratio": ratio, "status": status,
        })
    return pd.DataFrame(results)

# %% [markdown]
# ## 2. Liftoff gene table (apply the rename map)
#
# The rename list from `02_rename_loc_genes` gives non-`LOC` names to Liftoff
# genes. We apply it so the same-name ambiguity filter (step 6) and the paralog
# suffixing (step 7) see the final gene names.

# %%
gffl = pd.read_csv(LIFTOFF_GFF_GENES, sep="\t", comment="#", header=None,
                   names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
gffl_gene = gffl[gffl["type"] == "gene"].copy()
gffl_gene["gene"] = gffl_gene["attributes"].str.extract(r"ID=gene-([^;]+)", expand=False).str.upper()

rename_list = pd.read_csv(RENAME_LIST, sep="\t", index_col=0)
map_gene_list = dict(zip(
    rename_list["loc_gene"].str.removeprefix("gene-"),
    rename_list["new_gene_name"].str.replace(r"-l\d$", "", regex=True),
))
gffl_gene["gene_name"] = gffl_gene["gene"].map(map_gene_list).fillna(gffl_gene["gene"]).str.upper()

print(f"Liftoff genes: {len(gffl_gene)}  (still LOC: {(gffl_gene['gene_name'].str.startswith('LOC')).sum()})")

# %% [markdown]
# ## 3. Criterion 1 — no CDS overlap with Liftoff

# %%
gff_denovo = parse_sorted_gff(DENOVO_GFF)
gff_liftoff = parse_sorted_gff(LIFTOFF_GFF_CDS)

feature_noCDSOverlap = find_non_overlapping_genes_with_cds(gff_denovo, gff_liftoff)
print(f"De novo genes with no Liftoff CDS overlap: {feature_noCDSOverlap['gene_id'].nunique()}")

# %% [markdown]
# ## 4. Criterion 2 — CDS length ≥ 60 % of mouse homolog

# %%
mouse_cds = parse_mouse_gff(MOUSE_GFF)
completeness = filter_genes_with_exons(feature_noCDSOverlap, mouse_cds, MIN_MOUSE_CDS_RATIO)
print("Completeness summary:")
print(completeness["status"].value_counts().to_string())

complete_genes = completeness[completeness["status"] == "complete"]
feature_noCDSOverlap = feature_noCDSOverlap[feature_noCDSOverlap["gene_id"].isin(complete_genes["gene_id"])]
print(f"After mouse-CDS-length filter: {feature_noCDSOverlap['gene_id'].nunique()} genes")

# %% [markdown]
# ## 5. Criterion 3 — expression > 0

# %%
bw = pyBigWig.open(BIGWIG)
feature_noCDSOverlap = feature_noCDSOverlap.copy()
feature_noCDSOverlap["bw_mean"] = feature_noCDSOverlap.apply(
    lambda row: get_bw_stats(bw, row["chr"], row["start"], row["end"]), axis=1,
)
bw.close()

feature_noCDSOverlap = feature_noCDSOverlap[feature_noCDSOverlap["bw_mean"] > 0]
print(f"After expression filter (> 0): {feature_noCDSOverlap['gene_id'].nunique()} genes")

# %% [markdown]
# ## 6. Remove same-name gene-span ambiguities

# %%
feature_noCDSOverlap = feature_noCDSOverlap.copy()
feature_noCDSOverlap["gene_name_orig"] = feature_noCDSOverlap["gene_name"]
feature_noCDSOverlap["gene_name"] = feature_noCDSOverlap["gene_name"].str.upper()

n_before = len(feature_noCDSOverlap)
feature_noCDSOverlap = filter_entire_genes_with_overlaps(feature_noCDSOverlap, gffl_gene)
print(f"Filtered {n_before} -> {len(feature_noCDSOverlap)} CDS features (same-name gene-span overlap)")

# Drop de novo genes without a functional name.
feature_noCDSOverlap = feature_noCDSOverlap[~(feature_noCDSOverlap["gene_name"] == "")]

# %% [markdown]
# ## 7. Assign `-dl` suffix to paralogs

# %%
df = feature_noCDSOverlap.copy()
df["gene_name"] = df["gene_name"].str.split(",").str[0]

df["is_duplicated"] = df["gene_name"].duplicated(keep=False)
df["present_in_octDeg1"] = df["gene_name"].isin(gffl_gene["gene_name"])

df["new_gene_name"] = df["gene_name"]

# Case 1: duplicated and absent from OctDeg1 -> suffix copies 2, 3, ...
mask1 = (df["is_duplicated"] == True) & (df["present_in_octDeg1"] == False)
for gene in df.loc[mask1, "gene_name"].unique():
    indices = df[(df["gene_name"] == gene) & mask1].index
    for i, idx in enumerate(indices[1:], start=1):
        df.at[idx, "new_gene_name"] = f"{gene}-dl{i}"

# Case 2: duplicated and present in OctDeg1 -> suffix all copies.
mask2 = (df["is_duplicated"] == True) & (df["present_in_octDeg1"] == True)
for gene in df.loc[mask2, "gene_name"].unique():
    indices = df[(df["gene_name"] == gene) & mask2].index
    for i, idx in enumerate(indices, start=1):
        df.at[idx, "new_gene_name"] = f"{gene}-dl{i}"

# Case 3: single copy but present in OctDeg1 -> suffix the one copy.
mask3 = (df["is_duplicated"] == False) & (df["present_in_octDeg1"] == True)
for gene in df.loc[mask3, "gene_name"].unique():
    indices = df[(df["gene_name"] == gene) & mask3].index
    if len(indices) > 0:
        df.at[indices[0], "new_gene_name"] = f"{gene}-dl1"

# %% [markdown]
# ## 8. Write the add list

# %%
collapsed_df = df[["gene_id", "new_gene_name"]].drop_duplicates().reset_index(drop=True)
assert not collapsed_df["new_gene_name"].str.contains(",", na=False).any(), "names still contain commas"
collapsed_df.to_csv(OUT_ADD_LIST, sep="\t")
print(f"Wrote {len(collapsed_df)} genes to {OUT_ADD_LIST}")
print(collapsed_df.head())
