#!/usr/bin/env python3

import pandas as pd
import sys

# Input arguments
asv_mag_file = sys.argv[1]        # e.g. asv_to_mag_mapping_97.tsv
gtdbtk_file = sys.argv[2]         # e.g. all_samples_gtdbtk_summary.tsv
output_file = sys.argv[3]         # e.g. asv_to_mag_mapping_97_with_tax.tsv

# Load data
asv_df = pd.read_csv(asv_mag_file, sep="\t")
gtdb_df = pd.read_csv(gtdbtk_file, sep="\t")

# Clean GTDB-Tk bin name by removing 'Bin_' prefix
gtdb_df["bin"] = gtdb_df["user_genome"].str.replace(r"^Bin_", "", regex=True)

# Merge ASV-to-MAG mapping with GTDB taxonomy on both MAG_bin and MAG_sample
merged = asv_df.merge(
    gtdb_df,
    how="left",
    left_on=["MAG_bin", "MAG_sample"],
    right_on=["bin", "sample"]
)

# Drop helper 'bin' column used for merging
if "bin" in merged.columns:
    merged.drop(columns=["bin"], inplace=True)

# Remove duplicate rows caused by multiple GTDB matches (e.g. same MAG annotated multiple times)
merged.drop_duplicates(subset=["Cluster", "ASV", "MAG_bin", "MAG_sample"], inplace=True)

# Write final output
merged.to_csv(output_file, sep="\t", index=False)
