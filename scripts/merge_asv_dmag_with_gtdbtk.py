#!/usr/bin/env python3

import pandas as pd
import sys

# Input arguments
asv_dmag_file = sys.argv[1]        # e.g. asv_to_dmag_mapping.tsv
gtdbtk_file = sys.argv[2]          # e.g. dmag_gtdbtk_summary.tsv
output_file = sys.argv[3]          # e.g. asv_to_dmag_mapping_with_tax.tsv

# Load input files
asv_df = pd.read_csv(asv_dmag_file, sep="\t")
gtdb_df = pd.read_csv(gtdbtk_file, sep="\t")

# Merge ASV-dMAG mapping with GTDB-Tk using dMAG == user_genome
merged = asv_df.merge(
    gtdb_df,
    how="left",
    left_on="dMAG",
    right_on="user_genome"
)

# Drop redundant user_genome column (optional)
if "user_genome" in merged.columns:
    merged.drop(columns=["user_genome"], inplace=True)

# Drop potential duplicate rows
merged.drop_duplicates(subset=["Cluster", "ASV", "MAG_bin", "dMAG"], inplace=True)

# Write output
merged.to_csv(output_file, sep="\t", index=False)
