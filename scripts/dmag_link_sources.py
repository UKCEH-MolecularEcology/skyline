#!/usr/bin/env python3
"""
Per dereplicated MAG (dMAG): which MAGs carried the long-read 16S evidence behind its ASV/OTU links.

A link is ASV/OTU -> long-read 16S on a contig of MAG_bin -> MAG_bin's dRep winner (dMAG). MAG_bin can be the
winner itself or another member of its dRep cluster (which only an all-MAGs contig table can see).

Usage: dmag_link_sources.py <asv_to_dmag_mapping_with_tax.tsv> <Cdb.csv> <Wdb.csv> <out.tsv>

Output (one row per linked dMAG):
  n_linked_ASVs               ASVs/OTUs linked to the dMAG
  n_MAGs_with_16S_link        MAGs (winner and/or other cluster members) that carried a matching 16S
  winner_has_16S_link         True if the winner itself carried one
  n_other_MAGs_with_16S_link  other members of the dRep cluster that carried one
  drep_cluster_size           MAGs in the dRep cluster (Cdb)
  linked_MAGs                 the MAGs, ';'-separated
A one-line summary is printed (to the log).
"""
import sys

import pandas as pd

links_file, cdb_file, wdb_file, out_file = sys.argv[1:5]
strip_fa = lambda s: s.astype(str).str.strip().str.replace(r"\.(fa|fna|fasta)$", "", regex=True)

cdb = pd.read_csv(cdb_file)
wdb = pd.read_csv(wdb_file)
size = cdb.groupby("secondary_cluster").size()
drep_size = pd.Series(wdb["cluster"].map(size).values, index=strip_fa(wdb["genome"]))

cols = ["n_linked_ASVs", "n_MAGs_with_16S_link", "winner_has_16S_link",
        "n_other_MAGs_with_16S_link", "drep_cluster_size", "linked_MAGs"]
m = pd.read_csv(links_file, sep="\t", dtype=str)
if m.empty:
    pd.DataFrame(columns=["dMAG"] + cols).to_csv(out_file, sep="\t", index=False)
    print("0 dMAGs linked")
    sys.exit(0)
for c in ("ASV", "MAG_bin", "dMAG"):
    m[c] = m[c].str.strip()
m = m[m["dMAG"].notna() & (m["dMAG"] != "no_dmag")]

g = m.groupby("dMAG")
s = pd.DataFrame({
    "n_linked_ASVs": g["ASV"].nunique(),
    "n_MAGs_with_16S_link": g["MAG_bin"].nunique(),
    "winner_has_16S_link": (m["MAG_bin"] == m["dMAG"]).groupby(m["dMAG"]).any(),
})
s["n_other_MAGs_with_16S_link"] = s["n_MAGs_with_16S_link"] - s["winner_has_16S_link"].astype(int)
s["drep_cluster_size"] = drep_size.reindex(s.index).astype("Int64")
s["linked_MAGs"] = g["MAG_bin"].apply(lambda x: ";".join(sorted(set(x))))
s.index.name = "dMAG"
s.sort_values(["n_MAGs_with_16S_link", "n_linked_ASVs"], ascending=False)[cols].to_csv(out_file, sep="\t")

only_other = int((~s["winner_has_16S_link"]).sum())
print("{} dMAGs linked: {} via the winner's own 16S, {} only via other MAGs in their dRep cluster; "
      "MAGs with a 16S link per dMAG: median {:.0f}, max {}".format(
          len(s), int(s["winner_has_16S_link"].sum()), only_other,
          s["n_MAGs_with_16S_link"].median(), int(s["n_MAGs_with_16S_link"].max())))
