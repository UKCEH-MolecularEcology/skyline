#!/usr/bin/env python3

import pandas as pd
import sys
import re

# Input arguments
input_file = sys.argv[1]         # e.g., cluster_members.csv (comma-separated)
dmag_map_file = sys.argv[2]      # e.g., contig_MAGs_dMAGs.tsv (tab-separated)
output_file = sys.argv[3]        # e.g., asv_to_sampleMAG_dMAG.tsv

# Load cluster membership CSV
df = pd.read_csv(input_file)  # <-- assumes comma-separated

# Load dMAG mapping TSV
dmag_map = pd.read_csv(dmag_map_file, sep="\t")
valid_sample_mags = set(dmag_map["MAG"].unique())

# Extract sample-contig from long-read Members (e.g., b10_d3_1400_con-ctg11002996)
def extract_sample_contig(value):
    if pd.isnull(value):
        return None
    match = re.match(r"^(.*?-ctg\d+)", value)
    return match.group(1) if match else None

# Extract sample name from sample-contig (e.g., b10_d3_1400_con)
def extract_sample(sample_contig):
    if pd.isnull(sample_contig):
        return None
    return sample_contig.split("-ctg")[0]

# Add sample-contig and sample columns (only for long reads)
df["sample_contig"] = df.apply(
    lambda row: extract_sample_contig(row["Members"]) 
    if row["Member_Type"].lower() == "long_read" 
    else None,
    axis=1
)
df["sample"] = df["sample_contig"].apply(extract_sample)

# Filter valid MAGs that are linked to a dMAG
mag_info = df.loc[
    (df["Member_Type"].str.lower() == "long_read") &
    (df["MAG_member"].notnull()) &
    (df["MAG_member"] != "no_mag") &
    (df["MAG_member"].isin(valid_sample_mags)),
    ["Cluster", "MAG_member", "sample"]
].drop_duplicates().rename(columns={
    "MAG_member": "MAG_bin",
    "sample": "MAG_sample"
})

# Add dMAG name
mag_info = mag_info.merge(
    dmag_map[["MAG", "dMAG"]].drop_duplicates(),
    left_on="MAG_bin",
    right_on="MAG",
    how="left"
).drop(columns=["MAG"])

# Filter ASVs (amplicon-type members)
asv_info = df.loc[
    df["Member_Type"].str.lower().str.startswith("amplicon"),
    ["Cluster", "Members"]
].rename(columns={"Members": "ASV"})

# Join each ASV to each MAG in its cluster
results = []

for cluster_id, asvs in asv_info.groupby("Cluster"):
    mags = mag_info[mag_info["Cluster"] == cluster_id]
    for asv in asvs["ASV"]:
        for _, mag in mags.iterrows():
            results.append({
                "Cluster": cluster_id,
                "ASV": asv,
                "MAG_bin": mag["MAG_bin"],
                "MAG_sample": mag["MAG_sample"],
                "dMAG": mag["dMAG"]
            })

# Save result as tab-separated output
pd.DataFrame(results).to_csv(output_file, sep="\t", index=False)

#
#import pandas as pd
#import sys
#import re
#
## Input and output file paths
#input_file = sys.argv[1]
#output_file = sys.argv[2]
#
## Load input file
#df = pd.read_csv(input_file)
#
## Extract sample-contig from long-read Members (e.g., b10_d3_1400_con-ctg11002996)
#def extract_sample_contig(value):
#    if pd.isnull(value):
#        return None
#    match = re.match(r"^(.*?-ctg\d+)", value)
#    return match.group(1) if match else None
#
## Extract sample name from sample-contig (e.g., b10_d3_1400_con from b10_d3_1400_con-ctg11002996)
#def extract_sample(sample_contig):
#    if pd.isnull(sample_contig):
#        return None
#    return sample_contig.split("-ctg")[0]
#
## Add sample-contig and sample columns (only for long reads)
#df["sample_contig"] = df.apply(
#    lambda row: extract_sample_contig(row["Members"]) 
#    if row["Member_Type"].lower() == "long_read" 
#    else None,
#    axis=1
#)
#df["sample"] = df["sample_contig"].apply(extract_sample)
#
## Filter valid MAGs (long reads with valid MAGs, ignoring 'no_mag')
#mag_info = df.loc[
#    (df["Member_Type"].str.lower() == "long_read") &
#    (df["MAG_member"].notnull()) &
#    (df["MAG_member"] != "no_mag"),
#    ["Cluster", "MAG_member", "sample"]
#].drop_duplicates().rename(columns={
#    "MAG_member": "MAG_bin",
#    "sample": "MAG_sample"
#})
#
## Filter ASVs (amplicon-type members)
#asv_info = df.loc[
#    df["Member_Type"].str.lower().str.startswith("amplicon"),
#    ["Cluster", "Members"]
#].rename(columns={"Members": "ASV"})
#
## Join each ASV to each MAG in its cluster
#results = []
#
#for cluster_id, asvs in asv_info.groupby("Cluster"):
#    mags = mag_info[mag_info["Cluster"] == cluster_id]
#    for asv in asvs["ASV"]:
#        for _, mag in mags.iterrows():
#            results.append({
#                "Cluster": cluster_id,
#                "ASV": asv,
#                "MAG_bin": mag["MAG_bin"],
#                "MAG_sample": mag["MAG_sample"]
#            })
#
## Save result
#pd.DataFrame(results).to_csv(output_file, sep="\t", index=False)
