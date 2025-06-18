#!/usr/bin/env python3

import pandas as pd
import sys
import re

# Input arguments
bin_info_file = sys.argv[1]
cluster_file = sys.argv[2]
output_file = sys.argv[3]

# Load input files
bin_df = pd.read_csv(bin_info_file, sep="\t")
cluster_df = pd.read_csv(cluster_file)

# Correctly extract sample-contig (up to -ctgXXXXX) using regex
def extract_sample_contig(member):
    if pd.isnull(member):
        return None
    match = re.match(r"^(.*?-ctg\d+)", member)
    return match.group(1) if match else None

# Apply to extract normalized key for merging
cluster_df["normalized_member"] = cluster_df["Members"].apply(extract_sample_contig)

# Merge on normalized_member vs sample-contig
merged = cluster_df.merge(
    bin_df[["sample-contig", "bin"]],
    how="left",
    left_on="normalized_member",
    right_on="sample-contig"
)

# Add bin_matched column
merged["bin_matched"] = merged["bin"].notnull()

# Fill bin only for unmatched long reads
missing_bin_mask = (merged["bin"].isnull()) & (merged["Member_Type"].str.lower() == "long_read")
merged.loc[missing_bin_mask, "bin"] = "no_bin"

# Drop helper columns
# merged.drop(columns=["normalized_member", "sample-contig"], inplace=True)
merged.drop(columns=["sample-contig"], inplace=True)

# Write result
merged.to_csv(output_file, index=False)
