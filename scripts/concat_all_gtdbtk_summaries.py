#!/usr/bin/env python3

import pandas as pd
import sys
import os
from glob import glob

# Arguments
bin_info_dir = sys.argv[1]  # e.g. results/bin_info
output_file = sys.argv[2]

# Find all sample directories inside bin_info_dir
sample_dirs = [d for d in glob(os.path.join(bin_info_dir, "*")) if os.path.isdir(d)]

all_dfs = []

for sample_dir in sample_dirs:
    sample_name = os.path.basename(sample_dir)
    arc_file = os.path.join(sample_dir, f"{sample_name}_gtdbtk_arc_summary.tsv")
    bac_file = os.path.join(sample_dir, f"{sample_name}_gtdbtk_bac_summary.tsv")
    
    # Read and add sample column
    for file in [arc_file, bac_file]:
        if os.path.exists(file):
            df = pd.read_csv(file, sep="\t")
            df["sample"] = sample_name  # reinforce sample name from folder
            all_dfs.append(df)

# Concatenate all dataframes
combined_df = pd.concat(all_dfs, ignore_index=True)

# Move 'sample' column to front
cols = combined_df.columns.tolist()
cols.insert(0, cols.pop(cols.index('sample')))
combined_df = combined_df[cols]

# Save combined file
combined_df.to_csv(output_file, sep="\t", index=False)
