#!/usr/bin/env python3

import pandas as pd
import re
import sys

# Input files
mag_info_file = sys.argv[1]
cluster_file = sys.argv[2]
output_file = sys.argv[3]

# cluster_file = "results/OTU/99_cluster/clusters_99_mixed_clusters_split.csv"
# mag_info_file = "/ei/.project-scratch/5/542de014-1e71-4955-945a-5d2ab09567a7/CEHsoil/HiFi/assemblies/mags_with_rhyzo/results/contig_MAGs_dMAGs.tsv"
# output_file = "results/OTU/99_cluster/clusters_99_mixed_clusters_with_MAGs.tsv"

# Load data
cluster_df = pd.read_csv(cluster_file)
mag_df = pd.read_csv(mag_info_file, sep="\t")

# Extract normalized sample-contig ID (e.g. b10_d3_0200_con-ctg8980024)
def extract_sample_contig(val):
    if pd.isnull(val):
        return None
    m = re.match(r"^(.*?-ctg\d+)", val)
    return m.group(1) if m else None

# Normalize member and representative contig IDs
cluster_df["normalized_member"] = cluster_df["Members"].apply(extract_sample_contig)
cluster_df["normalized_representative"] = cluster_df["Representative"].apply(extract_sample_contig)

# Create sample-contig lookup in MAG dataframe
mag_df["sample-contig"] = mag_df["asm"] + "-" + mag_df["contig"]

# Merge MAG info on normalized member
merged = cluster_df.merge(
    mag_df[["sample-contig", "MAG", "dMAG"]],
    how="left",
    left_on="normalized_member",
    right_on="sample-contig"
)

merged.rename(columns={"MAG": "MAG_member", "dMAG": "dMAG_member"}, inplace=True)
merged.drop(columns=["sample-contig"], inplace=True)

# Fallback: merge using representative if MAG still missing
missing_mask = merged["MAG_member"].isnull()
if missing_mask.any():
    mag_lookup = mag_df.set_index("sample-contig")[["MAG", "dMAG"]]
    rep_mags = merged.loc[missing_mask, "normalized_representative"].map(mag_lookup["MAG"])
    rep_dmags = merged.loc[missing_mask, "normalized_representative"].map(mag_lookup["dMAG"])
    merged.loc[missing_mask, "MAG_member"] = rep_mags
    merged.loc[missing_mask, "dMAG_member"] = rep_dmags
    merged.loc[missing_mask, "mag_source"] = "Representative"

# Assign source if found through Member
merged.loc[merged["MAG_member"].notnull() & merged["mag_source"].isnull(), "mag_source"] = "Member"

# Fill in unmatched cases
merged["MAG_member"].fillna("no_mag", inplace=True)
merged["dMAG_member"].fillna("no_dmag", inplace=True)
merged["mag_source"].fillna("no_match", inplace=True)

# Output result
merged.to_csv(output_file, index=False)
