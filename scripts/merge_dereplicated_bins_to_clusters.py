#!/usr/bin/env python3
"""
Attach contig -> MAG -> dMAG information to each member of the mixed clusters.

A member gets a MAG only through its OWN contig (no representative fallback: that gave unbinned long-read 16S
the representative's MAG, with the member's assembly as provenance).

Usage: merge_dereplicated_bins_to_clusters.py <contig_MAGs_dMAGs.tsv> <clusters_<id>_mixed_clusters_split.csv> <out.csv>
Output columns: input columns + normalized_member, normalized_representative, MAG_member, dMAG_member, mag_source
"""
import re
import sys

import pandas as pd

mag_info_file, cluster_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

cluster_df = pd.read_csv(cluster_file)
mag_df = pd.read_csv(mag_info_file, sep="\t", dtype=str)


def extract_sample_contig(val):
    """<assembly>-ctg<digits> from a long-read 16S name (e.g. b10_d3_0200_con-ctg8980024_17408...18650)"""
    if pd.isnull(val):
        return None
    m = re.match(r"^(.*?-ctg\d+)", str(val).strip())
    return m.group(1) if m else None


cluster_df["normalized_member"] = cluster_df["Members"].apply(extract_sample_contig)
cluster_df["normalized_representative"] = cluster_df["Representative"].apply(extract_sample_contig)

mag_df["sample-contig"] = mag_df["asm"] + "-" + mag_df["contig"]
dup = mag_df["sample-contig"].duplicated(keep=False)
if dup.any():
    sys.exit("contig table: {} rows share an assembly-contig with another MAG, e.g. {}".format(
        int(dup.sum()), mag_df.loc[dup, "sample-contig"].iloc[0]))

merged = cluster_df.merge(mag_df[["sample-contig", "MAG", "dMAG"]], how="left",
                          left_on="normalized_member", right_on="sample-contig")
merged = merged.drop(columns=["sample-contig"]).rename(columns={"MAG": "MAG_member", "dMAG": "dMAG_member"})
merged["mag_source"] = merged["MAG_member"].notna().map({True: "Member", False: "no_match"}).astype(object)
merged["MAG_member"] = merged["MAG_member"].fillna("no_mag")
merged["dMAG_member"] = merged["dMAG_member"].fillna("no_dmag")
merged.to_csv(output_file, index=False)
print("{} members, {} with a MAG through their own contig".format(len(merged), int((merged["mag_source"] == "Member").sum())))
