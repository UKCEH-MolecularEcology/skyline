#!/usr/bin/env python3

import pandas as pd
import sys
import os

# Arguments
sample_dir = sys.argv[1]  # e.g. results/bin_info/r30_d3_1400_bio
output_file = sys.argv[2]

# Paths to the two GTDB-Tk summary files
arc_file = os.path.join(sample_dir, os.path.basename(sample_dir) + "_gtdbtk_arc_summary.tsv")
bac_file = os.path.join(sample_dir, os.path.basename(sample_dir) + "_gtdbtk_bac_summary.tsv")

# Read the files
arc_df = pd.read_csv(arc_file, sep="\t")
bac_df = pd.read_csv(bac_file, sep="\t")

# Concatenate
combined_df = pd.concat([arc_df, bac_df], ignore_index=True)

# Move 'sample' column to front
cols = combined_df.columns.tolist()
cols.insert(0, cols.pop(cols.index('sample')))
combined_df = combined_df[cols]

# Write output
combined_df.to_csv(output_file, sep="\t", index=False)
