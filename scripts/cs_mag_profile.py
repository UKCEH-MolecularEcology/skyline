#!/usr/bin/env python3
"""
Profile dereplicated MAGs across an amplicon survey (e.g. Countryside Survey 2007 16S) through 16S links.

Inputs
  --links      asv_to_dmag_mapping_<id>_with_tax.tsv (Cluster, ASV, MAG_bin, MAG_sample, dMAG, GTDB columns)
  --otu-table  samples x OTUs count table (CSV, first column = sample ID)
  --metadata   sample metadata (CSV), joined on --sample-column (e.g. PLOT)
  --taxonomy   optional OTU taxonomy (CSV: OTU id, taxonomy)
  --prefix     prefix that was added to the amplicon names for clustering (stripped to recover the OTU IDs)

Outputs (in --out-dir)
  otu_mag_links.tsv          one row per survey OTU: linked dMAG(s), GTDB + survey taxonomy, reads, prevalence
  link_summary.tsv           how many OTUs / reads link to MAGs and how many do not
  dmag_sample_abundance.tsv.gz   dMAG x sample: reads and relative abundance
                               reads_all    = all OTUs linked to the dMAG
                               reads_unique = only OTUs linked to this dMAG alone (no shared OTUs)
  dmag_habitat.tsv           dMAG x metadata category: prevalence and mean relative abundance
  dmag_soil.tsv              dMAG: mean of numeric metadata (e.g. pH) where present, and abundance-weighted
  dmag_summary.tsv           one row per linked dMAG: taxonomy, OTUs, samples present, top category per variable
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--links", required=True)
    ap.add_argument("--otu-table", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--taxonomy", default="")
    ap.add_argument("--prefix", default="")
    ap.add_argument("--sample-column", default="PLOT")
    ap.add_argument("--features-in-rows", action="store_true",
                    help="count table has features (ASVs/OTUs) as rows and samples as columns")
    ap.add_argument("--sample-name-replace", default="",
                    help="OLD:NEW applied to the count-table sample names, e.g. '.:-' (R turns '-' into '.')")
    ap.add_argument("--habitat-columns", nargs="*", default=[])
    ap.add_argument("--numeric-columns", nargs="*", default=[])
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    out = lambda name: os.path.join(a.out_dir, name)

    # ---- survey counts (samples x OTUs) and metadata
    counts = pd.read_csv(a.otu_table, index_col=0)
    if a.features_in_rows:
        counts = counts.T
    counts.index = counts.index.astype(str).str.strip()
    if a.sample_name_replace:
        old, new = a.sample_name_replace.split(":", 1)
        counts.index = counts.index.str.replace(old, new, regex=False)
    counts.columns = counts.columns.astype(str).str.strip()
    counts = counts.fillna(0)
    totals = counts.sum(axis=1)
    rel = counts.div(totals.replace(0, np.nan), axis=0).fillna(0)
    meta = pd.read_csv(a.metadata, dtype={a.sample_column: str})
    meta[a.sample_column] = meta[a.sample_column].str.strip()
    meta = meta.drop_duplicates(a.sample_column).set_index(a.sample_column)
    in_meta = counts.index.isin(meta.index)
    print("samples: {} in table, {} with metadata".format(len(counts), int(in_meta.sum())))

    # ---- links: amplicon name -> OTU id
    links = pd.read_csv(a.links, sep="\t", dtype=str)
    for c in ("ASV", "dMAG"):
        if c in links.columns:
            links[c] = links[c].str.strip()          # stray carriage returns (CRLF FASTA) / spaces in names
    if len(links):
        links["OTU"] = links["ASV"].str.replace("^" + a.prefix, "", regex=True) if a.prefix else links["ASV"]
        links = links[links["dMAG"].notna() & (links["dMAG"] != "no_dmag")]
    else:
        links = pd.DataFrame(columns=["OTU", "dMAG", "classification"])
    tax_col = "classification" if "classification" in links.columns else None
    pairs = links[["OTU", "dMAG"] + ([tax_col] if tax_col else [])].drop_duplicates(["OTU", "dMAG"])
    linked_not_in_table = pairs.loc[~pairs["OTU"].isin(counts.columns), "OTU"].nunique()
    pairs = pairs[pairs["OTU"].isin(counts.columns)]
    n_dmag_per_otu = pairs.groupby("OTU")["dMAG"].nunique()

    # ---- per-OTU table
    otu = pd.DataFrame(index=counts.columns)
    otu.index.name = "OTU"
    otu["reads"] = counts.sum(axis=0)
    otu["prevalence_samples"] = (counts > 0).sum(axis=0)
    otu["n_dMAGs"] = n_dmag_per_otu.reindex(otu.index).fillna(0).astype(int)
    otu["dMAGs"] = pairs.groupby("OTU")["dMAG"].apply(lambda x: ";".join(sorted(set(x)))).reindex(otu.index).fillna("")
    if tax_col:
        otu["dMAG_classification"] = pairs.groupby("OTU")[tax_col].apply(
            lambda x: ";".join(sorted(set(x.dropna())))).reindex(otu.index).fillna("")
    if a.taxonomy and os.path.exists(a.taxonomy):
        t = pd.read_csv(a.taxonomy, index_col=0)
        t.index = t.index.astype(str).str.strip()
        tax = t.iloc[:, 0] if t.shape[1] == 1 else t.fillna("").astype(str).agg("; ".join, axis=1)  # one column or ranks
        otu["survey_taxonomy"] = tax.reindex(otu.index)
    otu.sort_values(["n_dMAGs", "reads"], ascending=False).to_csv(out("otu_mag_links.tsv"), sep="\t")

    # ---- link summary
    linked = otu["n_dMAGs"] > 0
    lf = counts.loc[:, linked[linked].index].sum(axis=1) / totals.replace(0, np.nan)
    summ = [
        ("survey_OTUs", len(otu)),
        ("OTUs_linked_to_MAGs", int(linked.sum())),
        ("OTUs_not_linked", int((~linked).sum())),
        ("OTUs_linked_to_>1_dMAG", int((otu["n_dMAGs"] > 1).sum())),
        ("percent_OTUs_linked", round(100 * linked.mean(), 2)),
        ("percent_reads_linked_overall", round(100 * otu.loc[linked, "reads"].sum() / otu["reads"].sum(), 2)),
        ("percent_reads_linked_per_sample_median", round(100 * lf.median(), 2)),
        ("percent_reads_linked_per_sample_min", round(100 * lf.min(), 2)),
        ("percent_reads_linked_per_sample_max", round(100 * lf.max(), 2)),
        ("dMAGs_linked_to_survey_OTUs", int(pairs["dMAG"].nunique())),
        ("linked_OTUs_absent_from_count_table", int(linked_not_in_table)),
        ("samples", len(counts)),
        ("samples_with_metadata", int(in_meta.sum())),
    ]
    pd.DataFrame(summ, columns=["metric", "value"]).to_csv(out("link_summary.tsv"), sep="\t", index=False)
    for k, v in summ:
        print("  {:<42} {}".format(k, v))

    # ---- dMAG x sample abundance
    if pairs.empty:
        sys.stderr.write("no OTU-dMAG links: abundance/habitat tables will be empty\n")
    unique_otus = set(n_dmag_per_otu[n_dmag_per_otu == 1].index)
    rows = []
    for dmag, grp in pairs.groupby("dMAG"):
        otus = list(grp["OTU"])
        u = [o for o in otus if o in unique_otus]
        rows.append(pd.DataFrame({
            "dMAG": dmag, "sample": counts.index,
            "reads_all": counts[otus].sum(axis=1).values,
            "reads_unique": counts[u].sum(axis=1).values if u else 0,
            "rel_abund_all": rel[otus].sum(axis=1).values,
            "rel_abund_unique": rel[u].sum(axis=1).values if u else 0.0}))
    ab = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["dMAG", "sample", "reads_all", "reads_unique", "rel_abund_all", "rel_abund_unique"])
    ab[ab["reads_all"] > 0].to_csv(out("dmag_sample_abundance.tsv.gz"), sep="\t", index=False, compression="gzip")

    # ---- habitat / category profiles (all samples with metadata, including zeros)
    abm = ab.join(meta, on="sample", how="inner") if len(ab) else ab
    hab_rows = []
    for col in a.habitat_columns:
        if col not in meta.columns:
            sys.stderr.write("metadata column not found, skipped: {}\n".format(col))
            continue
        n_level = meta.loc[meta.index.isin(counts.index), col].value_counts()
        g = abm.assign(present=abm["reads_all"] > 0).groupby(["dMAG", col])
        h = g.agg(n_present=("present", "sum"), mean_rel_abund=("rel_abund_all", "mean")).reset_index()
        h["variable"] = col
        h = h.rename(columns={col: "level"})
        h["n_samples"] = h["level"].map(n_level)
        h["prevalence"] = h["n_present"] / h["n_samples"]
        hab_rows.append(h[["dMAG", "variable", "level", "n_samples", "n_present", "prevalence", "mean_rel_abund"]])
    hab = pd.concat(hab_rows, ignore_index=True) if hab_rows else pd.DataFrame(
        columns=["dMAG", "variable", "level", "n_samples", "n_present", "prevalence", "mean_rel_abund"])
    hab.sort_values(["dMAG", "variable", "prevalence"], ascending=[True, True, False]).to_csv(
        out("dmag_habitat.tsv"), sep="\t", index=False)

    # ---- numeric metadata (e.g. pH): where present, and abundance-weighted
    soil_rows = []
    for dmag, g in abm.groupby("dMAG") if len(abm) else []:
        pres = g[g["reads_all"] > 0]
        r = {"dMAG": dmag, "n_samples_present": len(pres)}
        for c in a.numeric_columns:
            if c not in g.columns:
                continue
            v = pd.to_numeric(pres[c], errors="coerce")
            w = pres["rel_abund_all"]
            ok = v.notna() & (w > 0)
            r[c + "_mean_present"] = v.mean()
            r[c + "_weighted"] = np.average(v[ok], weights=w[ok]) if ok.any() else np.nan
        soil_rows.append(r)
    pd.DataFrame(soil_rows).to_csv(out("dmag_soil.tsv"), sep="\t", index=False)

    # ---- per-dMAG summary
    if len(pairs):
        s = pairs.groupby("dMAG").agg(n_OTUs=("OTU", "nunique"))
        if tax_col:
            s["classification"] = pairs.groupby("dMAG")[tax_col].first()
        s["n_samples_present"] = ab[ab["reads_all"] > 0].groupby("dMAG")["sample"].nunique().reindex(s.index).fillna(0).astype(int)
        s["mean_rel_abund_all"] = ab.groupby("dMAG")["rel_abund_all"].mean().reindex(s.index)
        for col in a.habitat_columns:
            hc = hab[(hab["variable"] == col) & (hab["n_present"] > 0)]
            if len(hc):
                top = hc.sort_values("prevalence", ascending=False).groupby("dMAG").first()
                s["top_" + col] = top["level"].reindex(s.index)
                s["top_" + col + "_prevalence"] = top["prevalence"].reindex(s.index)
        s.sort_values("n_samples_present", ascending=False).to_csv(out("dmag_summary.tsv"), sep="\t")
    else:
        pd.DataFrame(columns=["dMAG"]).to_csv(out("dmag_summary.tsv"), sep="\t", index=False)
    print("done")


if __name__ == "__main__":
    main()
