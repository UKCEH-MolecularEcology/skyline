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

# Extract base contig ID using regex: keep only up to -ctgXXXXX
def extract_sample_contig(value):
    if pd.isnull(value):
        return None
    match = re.match(r"^(.*?-ctg\d+)", value)
    return match.group(1) if match else None

# Add normalized join keys for Members and Representative
cluster_df["normalized_member"] = cluster_df["Members"].apply(extract_sample_contig)
cluster_df["normalized_representative"] = cluster_df["Representative"].apply(extract_sample_contig)

# Merge on Members first
merged = cluster_df.merge(
    bin_df[["sample-contig", "bin"]],
    how="left",
    left_on="normalized_member",
    right_on="sample-contig"
)

# Initialize bin_source column
merged["bin_source"] = None
merged.loc[merged["bin"].notnull(), "bin_source"] = "Members"

# Fallback: fill in bin from Representative where still missing
missing_mask = merged["bin"].isnull()
if missing_mask.any():
    rep_lookup = bin_df.set_index("sample-contig")["bin"]
    merged.loc[missing_mask, "bin"] = merged.loc[missing_mask, "normalized_representative"].map(rep_lookup)
    # Update bin_source for those newly filled
    filled_from_rep = merged["bin_source"].isnull() & merged["bin"].notnull()
    merged.loc[filled_from_rep, "bin_source"] = "Representative"

# Mark whether we successfully matched a bin
merged["bin_matched"] = merged["bin"].notnull()

# Only fill no_bin for unmatched long reads
no_bin_mask = (merged["bin"].isnull()) & (merged["Member_Type"].str.lower() == "long_read")
merged.loc[no_bin_mask, "bin"] = "no_bin"
merged.loc[no_bin_mask, "bin_source"] = "no_bin"

# Drop helper columns
merged.drop(columns=["sample-contig"], inplace=True)

# Save to output
merged.to_csv(output_file, index=False)
