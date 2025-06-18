#!/usr/bin/env python3

import sys
import pandas as pd

# Inputs from Snakemake
circ_file = sys.argv[1]
consensus_file = sys.argv[2]
sample = sys.argv[3]
output_file = sys.argv[4]

# Read both clustering files
circ_df = pd.read_csv(circ_file)
cons_df = pd.read_csv(consensus_file)

# Concatenate them
combined = pd.concat([circ_df, cons_df], ignore_index=True)

# Rename column '0' to 'bin'
combined = combined.rename(columns={"0": "bin"})

# Add new column 'sample-contig'
combined["sample-contig"] = combined["sample"] + "-" + combined["contig_id"]

# Save to output
combined.to_csv(output_file, sep="\t", index=False)
