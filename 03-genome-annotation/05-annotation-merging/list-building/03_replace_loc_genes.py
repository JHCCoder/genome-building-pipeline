# %% [markdown]
# # Scenario 3 — Replace multiple Liftoff `LOC` genes with one de novo gene
#
# Replaces **multiple** Liftoff `LOC` genes with a **single** de novo (Braker3)
# gene when that de novo gene's CDS shows congruence (≥ 90 % same-strand overlap
# of each `LOC` CDS) with the CDS of several `LOC` genes — i.e. the de novo gene
# is a single, better model of a locus that Liftoff split into several `LOC`
# pieces.
#
# Paralogs are suffixed `-rl1`, `-rl2`, … so every gene name is unique.
#
# **Output:**
# * `gene_LOC_replace_unique.tsv` (`loc_gene<TAB>new_gene_name`)
# * `gene_LOC_replace_unique_with_annotated_gene.tsv` (`loc_gene<TAB>new_gene_name<TAB>annotated_gene`)
#
# **Run after** `02_rename_loc_genes` — its rename list is used to avoid name
# collisions when assigning `-rl` suffixes.

# %% [markdown]
# ## 0. Configuration

# %%
# --- inputs ---
LIFTOFF_GFF = "hifiasm-041425-scaffolded-chrAssigned-mito_agatProcessed.gff"  # Liftoff annotation (AGAT-processed)
DENOVO_GFF  = "braker_3UTRincluded_agatProcessed.gff"                         # Braker3 de novo (functionally named)
RENAME_LIST = "gene_LOC_nameChange_unique.tsv"                                # produced by 02_rename_loc_genes
MIN_CDS_OVERLAP = 0.9   # congruence threshold: fraction of the LOC CDS covered

# --- outputs ---
OUT_REPLACE_LIST        = "gene_LOC_replace_unique.tsv"
OUT_REPLACE_LIST_GENE   = "gene_LOC_replace_unique_with_annotated_gene.tsv"

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


def set_replace_new_gene_name(dfr, df=None):
    """Assign `-rl` suffixes to replaced genes so names stay unique.

    ``df`` (optional) is the rename list, used to avoid colliding with names
    already assigned by scenario 2.
    """
    dfr = dfr.copy()
    dfr["new_gene_name"] = dfr["annotated_mrna_upper"]

    # Case 1: name already present in OctDeg1 -> suffix every copy.
    for name, group in dfr[dfr["present_in_octDeg1"]].groupby("annotated_mrna_upper"):
        if any(group["loc_upper"].str.contains(name, regex=False)):
            dfr.loc[group.index, "new_gene_name"] = name
        else:
            for i, idx in enumerate(group.index, start=1):
                dfr.loc[idx, "new_gene_name"] = f"{name}-rl{i}"

    # Case 2: name absent from OctDeg1 -> suffix only if it collides with a renamed gene.
    if df is not None:
        for idx, row in dfr[~dfr["present_in_octDeg1"]].iterrows():
            name = row["annotated_mrna_upper"]
            if df["new_gene_name"].str.contains(name, regex=False).any():
                dfr.loc[idx, "new_gene_name"] = f"{name}-rl1"
            else:
                dfr.loc[idx, "new_gene_name"] = name
    else:
        mask = ~dfr["present_in_octDeg1"]
        dfr.loc[mask, "new_gene_name"] = dfr.loc[mask, "annotated_mrna_upper"]

    return dfr

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

cds_overlaps = find_cds_overlaps(df_liftoff, df_denovo, min_overlap=MIN_CDS_OVERLAP)
print(f"CDS pairs with >= {MIN_CDS_OVERLAP:.0%} same-strand overlap: {len(cds_overlaps)}")

# %% [markdown]
# ## 3. Annotate the overlap table

# %%
cds_overlaps["loc_mrna_name_upper"] = cds_overlaps["loc_gene"].str.extract(r"gene-([^;]+)", expand=False).str.upper()

mask = cds_overlaps["ref_mrna_name"].isna() | (cds_overlaps["ref_mrna_name"] == "")
cds_overlaps.loc[mask, "ref_mrna_name"] = cds_overlaps.loc[mask, "ref_cds_id"].str.split(".").str[0]

cds_overlaps["ref_mrna_name_upper"] = cds_overlaps["ref_mrna_name"].str.upper()
cds_overlaps["present_in_octDeg1"] = cds_overlaps["ref_mrna_name_upper"].isin(gffl_gene["gene"])

filtered_df1 = cds_overlaps[cds_overlaps["loc_gene"].str.fullmatch(r"gene-LOC\d{9}")].copy()
filtered_df = cds_overlaps[cds_overlaps["annotated_gene"].isin(filtered_df1["annotated_gene"])].copy()

filtered_df["num_annotated_genes_overlapping_loc_gene"] = filtered_df.groupby("loc_gene")["annotated_gene"].transform("nunique")
filtered_df["num_loc_genes_overlapping_annotated_gene"] = filtered_df.groupby("annotated_gene")["loc_gene"].transform("nunique")

# %% [markdown]
# ## 4. Replace decision — one de novo gene congruent with multiple `LOC` genes

# %%
rows_to_replace = []
for _, group in filtered_df.groupby("annotated_gene"):
    # The de novo gene must be annotated in every row.
    if not (group["annotated_gene"] != group["ref_mrna_name"]).all():
        continue
    multiple_loc = group["loc_gene"].nunique() > 1                      # de novo -> >1 LOC
    loc_or_same = group.apply(                                          # each LOC is LOC... or already the same name
        lambda row: bool(re.fullmatch(r"LOC\d{9}", row["loc_mrna_name_upper"]))
                    or (row["loc_mrna_name_upper"] == row["ref_mrna_name_upper"]),
        axis=1,
    ).all()
    many_to_one = (group["num_loc_genes_overlapping_annotated_gene"] > 1).all()  # de novo -> >1 LOC
    loc_to_one = (group["num_annotated_genes_overlapping_loc_gene"] == 1).all()  # each LOC -> this de novo only
    if multiple_loc and loc_or_same and many_to_one and loc_to_one:
        rows_to_replace.append(group)

replace_genes = pd.concat(rows_to_replace) if rows_to_replace else pd.DataFrame()
filtered_df["rename_replace"] = "no"
if not replace_genes.empty:
    filtered_df.loc[filtered_df["annotated_gene"].isin(replace_genes["annotated_gene"]), "rename_replace"] = "replace"
print(f"De novo genes replacing multiple LOC genes: {filtered_df[filtered_df['rename_replace'] == 'replace']['annotated_gene'].nunique()}")

# %% [markdown]
# ## 5. Build the replace table and assign `-rl` suffixes

# %%
change_df = filtered_df[filtered_df["rename_replace"] == "replace"].copy()
change_df = change_df[change_df["annotated_gene"] != change_df["ref_mrna_name"]]

replace_df = change_df[["loc_gene", "ref_mrna_name_upper", "annotated_gene",
                        "present_in_octDeg1", "loc_mrna_name_upper"]].drop_duplicates().reset_index(drop=True)
replace_df = replace_df.rename(columns={"ref_mrna_name_upper": "annotated_mrna_upper"})
replace_df["loc_upper"] = replace_df["loc_gene"].str.extract(r"gene-([^;]+)", expand=False).str.upper()

# Drop de novo genes that resolve to more than one annotated name (inconsistent).
inconsistent = replace_df.groupby("annotated_gene")["annotated_mrna_upper"].nunique()
inconsistent = inconsistent[inconsistent > 1].index.tolist()
if inconsistent:
    replace_df = replace_df[~replace_df["annotated_gene"].isin(inconsistent)].reset_index(drop=True)

rename_list = pd.read_csv(RENAME_LIST, sep="\t", index_col=0)  # columns: loc_gene, new_gene_name
replace_df = set_replace_new_gene_name(replace_df, df=rename_list)

# %% [markdown]
# ## 6. Write the replace lists

# %%
replace_df[["loc_gene", "new_gene_name"]].to_csv(OUT_REPLACE_LIST, sep="\t")
replace_df[["loc_gene", "new_gene_name", "annotated_gene"]].to_csv(OUT_REPLACE_LIST_GENE, sep="\t")
print(f"Wrote {len(replace_df)} genes to {OUT_REPLACE_LIST}")
print(replace_df[["loc_gene", "new_gene_name", "annotated_gene"]].head())
