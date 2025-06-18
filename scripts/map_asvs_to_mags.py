#!/usr/bin/env python3

import pandas as pd
import sys
import re

input_file = sys.argv[1]
output_file = sys.argv[2]

# Load input file
df = pd.read_csv(input_file)

# Extract sample-contig from MAG entries
def extract_sample_contig(value):
    if pd.isnull(value):
        return None
    match = re.match(r"^(.*?-ctg\d+)", value)
    return match.group(1) if match else None

# Extract sample name from sample-contig (before -ctg)
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

# Filter valid MAGs (long reads with bin info)
mag_info = df.loc[
    (df["Member_Type"].str.lower() == "long_read") &
    (df["bin"].notnull()) &
    (df["bin"] != "no_bin"),
    ["Cluster", "bin", "sample"]
].drop_duplicates().rename(columns={"bin": "MAG_bin", "sample": "MAG_sample"})

# Filter ASVs
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
                "MAG_sample": mag["MAG_sample"]
            })

# Save result
pd.DataFrame(results).to_csv(output_file, sep="\t", index=False)
