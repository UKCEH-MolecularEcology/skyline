#!/usr/bin/env python3

import pandas as pd
import sys

# Input files
asv_mag_file = sys.argv[1]        # e.g. asv_to_mag_mapping_97.tsv
gtdbtk_file = sys.argv[2]         # e.g. all_samples_gtdbtk_summary.tsv
output_file = sys.argv[3]         # e.g. asv_to_mag_mapping_97_with_tax.tsv

# Load data
asv_df = pd.read_csv(asv_mag_file, sep="\t")
gtdb_df = pd.read_csv(gtdbtk_file, sep="\t")

# Clean GTDB-Tk bin name
gtdb_df["bin"] = gtdb_df["user_genome"].str.replace(r"^Bin_", "", regex=True)

# Merge on MAG_bin <-> bin and MAG_sample <-> sample
merged = asv_df.merge(
    gtdb_df,
    how="left",
    left_on=["MAG_bin", "MAG_sample"],
    right_on=["bin", "sample"]
)

# Drop helper 'bin' column if redundant
merged.drop(columns=["bin"], inplace=True)

# Write merged output
merged.to_csv(output_file, sep="\t", index=False)
