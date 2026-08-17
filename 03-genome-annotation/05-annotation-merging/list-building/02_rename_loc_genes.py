# %% [markdown]
# # Scenario 2 — Rename Liftoff `LOC` genes
#
# Renames an unannotated (`LOC##########`) Liftoff gene to a de novo (Braker3)
# gene name when the two are the **same locus**. Congruence is defined as a
# **≥ 90 % overlap of the Liftoff CDS by the de novo CDS on the same strand**
# (one-directional overlap of the LOC CDS).
#
# A `LOC` gene is renamed when its CDS shows **unique** congruence — i.e. all of
# its overlapping, annotated de novo partners resolve to a *single* de novo gene
# name, and the mapping is 1:1 in both directions (the de novo gene overlaps only
# this one `LOC`, and this `LOC` overlaps only that one de novo gene).
#
# Paralogs are suffixed `-l1`, `-l2`, … so every gene name is unique.
#
# **Output:** `gene_LOC_nameChange_unique.tsv` (`loc_gene<TAB>new_gene_name`).
#
# **Run before** `01_add_novel_genes` (which consumes this rename list) and
# `03_replace_loc_genes` (which uses it to avoid name collisions).

# %% [markdown]
# ## 0. Configuration

# %%
# --- inputs ---
LIFTOFF_GFF = "hifiasm-041425-scaffolded-chrAssigned-mito_agatProcessed.gff"  # Liftoff annotation (AGAT-processed)
DENOVO_GFF  = "braker_3UTRincluded_agatProcessed.gff"                         # Braker3 de novo (functionally named)
MIN_CDS_OVERLAP = 0.9   # congruence threshold: fraction of the LOC CDS covered

# --- outputs ---
OUT_RENAME_LIST = "gene_LOC_nameChange_unique.tsv"

# %% [markdown]
# ## 1. Imports and helper functions

# %%
import pandas as pd
import numpy as np
import re


def parse_sorted_gff(gff_file):
    """Fast parser for a *sorted* GFF3 file (gene -> mRNA -> CDS/exon)."""
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
                    "strand": strand, "phase": phase, "attributes": attr_dict, "children": {},
                }
                current_mrna = None
            elif feature_type == "mRNA":
                if current_gene is None:
                    continue
                current_mrna = attr_dict["ID"]
                gene_dict[current_gene]["children"][current_mrna] = {
                    "type": "mRNA", "seqid": seqid, "source": source,
                    "start": int(start), "end": int(end), "score": score,
                    "strand": strand, "phase": phase, "attributes": attr_dict, "children": {},
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
                if feature_type == "CDS":
                    feature_id = f"{attr_dict['ID']}|{start}-{end}"
                gene_dict[current_gene]["children"][current_mrna]["children"][feature_id] = feature_data
    return gene_dict


def extract_cds(gff_dict, source_gff):
    """Flatten CDS features of a parsed GFF into a DataFrame."""
    rows = []
    for gene_id, gene_data in gff_dict.items():
        for mrna_id, mrna_data in gene_data["children"].items():
            mrna_attrs = mrna_data.get("attributes", {})
            protein_name = mrna_attrs.get("Name", "")
            for feature_id, feature_data in mrna_data["children"].items():
                if feature_data["type"] == "CDS":
                    rows.append({
                        "gene_id": gene_id,
                        "gene_name": gene_data["attributes"].get("Name", ""),
                        "mrna_id": mrna_id,
                        "mrna_name": protein_name,
                        "cds_id": feature_id,
                        "chr": feature_data["seqid"],
                        "start": int(feature_data["start"]),
                        "end": int(feature_data["end"]),
                        "strand": feature_data["strand"],
                        "phase": feature_data.get("phase", "."),
                        "source_gff": source_gff,
                    })
    return pd.DataFrame(rows)


def find_cds_overlaps(df1, df2, min_overlap=0.8):
    """Find CDS pairs on the same chromosome/strand where df1's CDS is ≥
    ``min_overlap`` covered by df2's CDS.

    df1 = Liftoff (`LOC`) CDS, df2 = de novo (annotated) CDS.
    """
    overlaps = []
    groups1 = df1.groupby(["chr", "strand"])
    groups2 = df2.groupby(["chr", "strand"])
    for (chr_val, strand_val), group1 in groups1:
        if (chr_val, strand_val) not in groups2.groups:
            continue
        group2 = groups2.get_group((chr_val, strand_val))
        overlap_starts = np.maximum.outer(group1["start"].values, group2["start"].values)
        overlap_ends = np.minimum.outer(group1["end"].values, group2["end"].values)
        overlap_lengths = overlap_ends - overlap_starts + 1
        overlap_lengths[overlap_lengths < 0] = 0
        lengths1 = (group1["end"].values - group1["start"].values + 1)[:, None]
        pct_overlaps = overlap_lengths / lengths1
        rows, cols = (pct_overlaps >= min_overlap).nonzero()
        for r, c in zip(rows, cols):
            cds1, cds2 = group1.iloc[r], group2.iloc[c]
            overlaps.append({
                "loc_gene": cds1["gene_id"],
                "annotated_gene": cds2["gene_id"],
                "loc_cds_id": cds1["cds_id"],
                "ref_cds_id": cds2["cds_id"],
                "ref_mrna_name": cds2.get("mrna_name", cds2.get("gene_name", "")),
                "loc_chr": cds1["chr"], "loc_start": cds1["start"], "loc_end": cds1["end"],
                "ref_start": cds2["start"], "ref_end": cds2["end"],
                "strand": cds1["strand"],
                "pct_overlap": pct_overlaps[r, c],
            })
    return pd.DataFrame(overlaps)

# %% [markdown]
# ## 2. Load annotations and find congruent CDS pairs

# %%
gffl = pd.read_csv(LIFTOFF_GFF, sep="\t", comment="#", header=None,
                   names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
gffl_gene = gffl[gffl["type"] == "gene"].copy()
gffl_gene["gene"] = gffl_gene["attributes"].str.extract(r"ID=gene-([^;]+)", expand=False).str.upper()

gff_denovo = parse_sorted_gff(DENOVO_GFF)
gff_liftoff = parse_sorted_gff(LIFTOFF_GFF)

df_denovo = extract_cds(gff_denovo, "gff_denovo")
df_liftoff = extract_cds(gff_liftoff, "gff_liftOff")

# df1 = Liftoff (LOC), df2 = de novo (annotated).
cds_overlaps = find_cds_overlaps(df_liftoff, df_denovo, min_overlap=MIN_CDS_OVERLAP)
print(f"CDS pairs with >= {MIN_CDS_OVERLAP:.0%} same-strand overlap: {len(cds_overlaps)}")

# %% [markdown]
# ## 3. Annotate the overlap table

# %%
# LOC gene name (upper), e.g. gene-LOC111814891 -> LOC111814891.
cds_overlaps["loc_mrna_name_upper"] = cds_overlaps["loc_gene"].str.extract(r"gene-([^;]+)", expand=False).str.upper()

# Fill empty de novo names from the de novo CDS id (gene id part).
mask = cds_overlaps["ref_mrna_name"].isna() | (cds_overlaps["ref_mrna_name"] == "")
cds_overlaps.loc[mask, "ref_mrna_name"] = cds_overlaps.loc[mask, "ref_cds_id"].str.split(".").str[0]

cds_overlaps["ref_mrna_name_upper"] = cds_overlaps["ref_mrna_name"].str.upper()
cds_overlaps["present_in_octDeg1"] = cds_overlaps["ref_mrna_name_upper"].isin(gffl_gene["gene"])

# Keep LOC genes and their overlapping de novo partners.
filtered_df1 = cds_overlaps[cds_overlaps["loc_gene"].str.fullmatch(r"gene-LOC\d{9}")].copy()
filtered_df = cds_overlaps[cds_overlaps["annotated_gene"].isin(filtered_df1["annotated_gene"])].copy()

# How many de novo genes does each LOC overlap, and vice versa.
filtered_df["num_annotated_genes_overlapping_loc_gene"] = filtered_df.groupby("loc_gene")["annotated_gene"].transform("nunique")
filtered_df["num_loc_genes_overlapping_annotated_gene"] = filtered_df.groupby("annotated_gene")["loc_gene"].transform("nunique")

# %% [markdown]
# ## 4. Rename decision — unique (1:1) congruence

# %%
rows_to_rename = []
for _, group in filtered_df.groupby("loc_gene"):
    # Only consider de novo partners that are actually annotated (have a name).
    g = group[group["annotated_gene"] != group["ref_mrna_name"]]
    if len(g) == 0:
        continue
    one_name = g["ref_mrna_name"].nunique() == 1                       # one unique de novo name
    one_to_one = (g["num_loc_genes_overlapping_annotated_gene"] == 1).all()   # de novo -> 1 LOC
    loc_to_one = (g["num_annotated_genes_overlapping_loc_gene"] == 1).all()   # LOC -> 1 de novo
    if one_name and one_to_one and loc_to_one:
        rows_to_rename.append(g)

rename_genes = pd.concat(rows_to_rename) if rows_to_rename else pd.DataFrame()
filtered_df["rename_replace"] = "no"
filtered_df.loc[filtered_df["loc_gene"].isin(rename_genes["loc_gene"]), "rename_replace"] = "rename"
print(f"LOC genes to rename: {filtered_df[filtered_df['rename_replace'] == 'rename']['loc_gene'].nunique()}")

# %% [markdown]
# ## 5. Assign `-l` suffix to paralogs

# %%
change_df = filtered_df[filtered_df["rename_replace"] == "rename"].copy()
# Drop rows whose de novo partner is unannotated (annotated_gene == ref_mrna_name).
change_df = change_df[change_df["annotated_gene"] != change_df["ref_mrna_name"]]

rename_df = change_df[["loc_gene", "ref_mrna_name_upper", "present_in_octDeg1"]].drop_duplicates().reset_index(drop=True)
rename_df["is_duplicated"] = rename_df["ref_mrna_name_upper"].duplicated(keep=False)

df = rename_df.copy()
df["new_gene_name"] = df["ref_mrna_name_upper"]

# Case 1: duplicated and absent from OctDeg1 -> suffix copies 2, 3, ...
mask1 = (df["is_duplicated"] == True) & (df["present_in_octDeg1"] == False)
for gene in df.loc[mask1, "ref_mrna_name_upper"].unique():
    indices = df[(df["ref_mrna_name_upper"] == gene) & mask1].index
    for i, idx in enumerate(indices[1:], start=1):
        df.at[idx, "new_gene_name"] = f"{gene}-l{i}"

# Case 2: present in OctDeg1 -> suffix all copies (name would collide with the reference).
mask2 = df["present_in_octDeg1"] == True
for gene in df.loc[mask2, "ref_mrna_name_upper"].unique():
    indices = df[(df["ref_mrna_name_upper"] == gene) & mask2].index
    for i, idx in enumerate(indices, start=1):
        df.at[idx, "new_gene_name"] = f"{gene}-l{i}"

# %% [markdown]
# ## 6. Write the rename list

# %%
rename_out = df[["loc_gene", "new_gene_name"]]
rename_out.to_csv(OUT_RENAME_LIST, sep="\t")  # default index column matches the merge tool / add notebook
print(f"Wrote {len(rename_out)} genes to {OUT_RENAME_LIST}")
print(rename_out.head())
