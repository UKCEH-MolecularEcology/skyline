#!/usr/bin/env python3
"""
Merge all MAG annotations into one table: one row per MAG, column groups by prefix.

  coverm_*        CoverM genome summary (prevalence, mean/max relative abundance, reads)
  abund_<sample>  relative abundance per sample (CoverM matrix; all samples kept)
  mt_*            microTrait traits (granularity --mt-granularity) + mt_mingentime, mt_optimumT
  pt_<tool>|<trait>[|<model>]_p / _b   porTraits probability / binary; pt_recognise_specI
  rgi_*           RGI: hits, Perfect, Strict, unique AROs, drug classes, gene families
  bacmet_*        BacMet: hits, genes
  as_*            antiSMASH: regions, complete (not on contig edge), with GECCO overlap, per class
  gc_*            GECCO: clusters, per type
  bs_*            BiG-SLiCE: GCFs, GCFs found only in this MAG, GCF ids

Usage:
  python merge_annotations.py --results-dir <annotations> --mags-dir <dereplicated_genomes> [--mags-ext fa]
         [--out <annotations>/merged/skyline_MAG_annotations.tsv.gz]
Missing inputs are skipped with a warning; counts are 0 for MAGs without hits.
"""
import argparse, glob, gzip, os, sys
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--results-dir", required=True)
p.add_argument("--mags-dir", required=True)
p.add_argument("--mags-ext", default="fa")
p.add_argument("--abund-matrix", default="coverm/coverm_rel_abund_masked_cf0.1.tsv",
               help="relative to --results-dir")
p.add_argument("--mt-granularity", default="3", choices=["1", "2", "3"])
p.add_argument("--out", default=None)
a = p.parse_args()
R = a.results_dir
out = a.out or os.path.join(R, "merged", "skyline_MAG_annotations.tsv.gz")
os.makedirs(os.path.dirname(out), exist_ok=True)

def path(rel):
    f = os.path.join(R, rel)
    if not os.path.exists(f):
        print("WARNING: missing, skipped:", f, file=sys.stderr); return None
    return f

def uniq_join(s):
    return ";".join(sorted({x for v in s.dropna().astype(str) for x in v.split(";") if x and x != "nan"}))

def safe(name):
    return str(name).strip().replace("\t", " ").replace(" ", "_")

mags = sorted(os.path.basename(f)[:-(len(a.mags_ext) + 1)] for f in glob.glob(os.path.join(a.mags_dir, "*." + a.mags_ext)))
M = pd.DataFrame(index=pd.Index(mags, name="Genome"))
print("MAGs:", len(M))
blocks, counts_cols = [], []

# ---- CoverM
f = path("coverm/coverm_genome_summary.tsv")
if f:
    d = pd.read_csv(f, sep="\t", index_col=0)
    blocks.append(d.add_prefix("coverm_"))
f = path(a.abund_matrix)
if f:
    d = pd.read_csv(f, sep="\t", index_col=0)
    blocks.append(d.add_prefix("abund_"))
    print("samples in abundance matrix:", d.shape[1])

# ---- microTrait
f = path("microtrait/tables/trait_matrixatgranularity{}.tsv".format(a.mt_granularity))
if f:
    d = pd.read_csv(f, sep="\t").set_index("id")
    extra = [c for c in ("mingentime", "optimumT") if c in d.columns]
    d = d[[c for c in d.columns if c not in extra] + extra]
    d.columns = ["mt_" + safe(c) for c in d.columns]
    blocks.append(d)

# ---- porTraits (wide: probability + binary per tool x trait [x model if several])
f = path("portraits/portraits_results.tsv.gz")
if f:
    d = pd.read_csv(f, sep="\t", dtype={"value_probability": float, "value_binary": float})
    n_models = d.groupby(["tool", "feature"])["tool_feature"].transform("nunique")
    d["key"] = "pt_" + d["tool"].astype(str) + "|" + d["feature"].astype(str)
    d.loc[n_models > 1, "key"] = d["key"] + "|" + d["tool_feature"].astype(str)
    prob = d.pivot_table(index="genome", columns="key", values="value_probability", aggfunc="first")
    binr = d.pivot_table(index="genome", columns="key", values="value_binary", aggfunc="first")
    pt = pd.concat([prob.add_suffix("_p"), binr.add_suffix("_b")], axis=1)
    pt = pt[sorted(pt.columns)]
    spec = d.groupby("genome")["speci"].first().rename("pt_recognise_specI")
    blocks += [pt, spec.to_frame()]
    print("porTraits columns:", pt.shape[1])

# ---- RGI
f = path("amr/rgi_all.tsv")
if f:
    d = pd.read_csv(f, sep="\t", low_memory=False)
    g = d.groupby("Genome")
    rgi = pd.DataFrame({
        "rgi_n_hits": g.size(),
        "rgi_n_perfect": g["Cut_Off"].apply(lambda s: (s == "Perfect").sum()),
        "rgi_n_strict": g["Cut_Off"].apply(lambda s: (s == "Strict").sum()),
        "rgi_n_unique_aro": g["Best_Hit_ARO"].nunique(),
        "rgi_genes": g["Best_Hit_ARO"].apply(uniq_join),
        "rgi_drug_classes": g["Drug Class"].apply(uniq_join) if "Drug Class" in d else None,
        "rgi_gene_families": g["AMR Gene Family"].apply(uniq_join) if "AMR Gene Family" in d else None,
    })
    blocks.append(rgi); counts_cols += ["rgi_n_hits", "rgi_n_perfect", "rgi_n_strict", "rgi_n_unique_aro"]

# ---- BacMet
f = path("amr/bacmet_all.tsv")
if f:
    d = pd.read_csv(f, sep="\t")
    g = d.groupby("Genome")
    bm = pd.DataFrame({"bacmet_n_hits": g.size(), "bacmet_n_unique_genes": g["bacmet_gene"].nunique(),
                       "bacmet_genes": g["bacmet_gene"].apply(uniq_join)})
    blocks.append(bm); counts_cols += ["bacmet_n_hits", "bacmet_n_unique_genes"]

# ---- antiSMASH (+ GECCO confirmation from the overlap table)
f = path("bgc/antismash_regions_all.tsv")
if f:
    d = pd.read_csv(f, sep="\t")
    d["complete"] = ~d["contig_edge"].astype(str).str.lower().eq("true")
    g = d.groupby("Genome")
    asm = pd.DataFrame({"as_n_regions": g.size(), "as_n_complete": g["complete"].sum(),
                        "as_products": g["products"].apply(uniq_join)})
    cls = d.assign(product=d["products"].astype(str).str.split(";")).explode("product").reset_index(drop=True)
    cls = pd.crosstab(cls["Genome"], cls["product"]).add_prefix("as_class_")
    cls.columns = [safe(c) for c in cls.columns]
    blocks += [asm, cls]; counts_cols += ["as_n_regions", "as_n_complete"] + list(cls.columns)
    fo = path("bgc/bgc_overlap.tsv")
    if fo:
        o = pd.read_csv(fo, sep="\t", usecols=["status", "Genome", "antismash_region"])
        conf = o[o["status"] == "both"].groupby("Genome")["antismash_region"].nunique().rename("as_n_with_gecco")
        blocks.append(conf.to_frame()); counts_cols.append("as_n_with_gecco")

# ---- GECCO
f = path("bgc/gecco_clusters_all.tsv")
if f:
    d = pd.read_csv(f, sep="\t")
    gc = pd.DataFrame({"gc_n_clusters": d.groupby("Genome").size()})
    types = d.assign(t=d["type"].astype(str).str.split(";")).explode("t").reset_index(drop=True)
    types = pd.crosstab(types["Genome"], types["t"]).add_prefix("gc_type_")
    types.columns = [safe(c) for c in types.columns]
    blocks += [gc, types]; counts_cols += ["gc_n_clusters"] + list(types.columns)

# ---- BiG-SLiCE
f = path("bgc/bigslice/tables/gcf_membership.tsv")
if f:
    d = pd.read_csv(f, sep="\t", usecols=["Genome", "gcf"])
    n_mag_per_gcf = d.groupby("gcf")["Genome"].nunique()
    d["singleton"] = d["gcf"].map(n_mag_per_gcf).eq(1)
    g = d.groupby("Genome")
    bs = pd.DataFrame({"bs_n_gcfs": g["gcf"].nunique(),
                       "bs_n_gcfs_only_this_mag": g.apply(lambda x: x.loc[x["singleton"], "gcf"].nunique()),
                       "bs_gcfs": g["gcf"].apply(lambda s: ";".join(sorted(set(s))))})
    blocks.append(bs); counts_cols += ["bs_n_gcfs", "bs_n_gcfs_only_this_mag"]

# ---- merge
for b in blocks:
    b.index = b.index.astype(str)
    missing = set(b.index) - set(M.index)
    if missing:
        print("WARNING: {} ids not among the MAGs (e.g. {}) in block starting {}".format(
            len(missing), next(iter(missing)), b.columns[0]), file=sys.stderr)
M = M.join(blocks, how="left") if blocks else M
for c in counts_cols:
    if c in M:
        M[c] = M[c].fillna(0).astype(int)
M.to_csv(out, sep="\t", compression="gzip")
print("wrote {}: {} MAGs x {} columns".format(out, M.shape[0], M.shape[1]))

# column dictionary
groups = [("coverm_", "CoverM genome summary"), ("abund_", "relative abundance in sample (" + a.abund_matrix + ")"),
          ("mt_", "microTrait granularity " + a.mt_granularity + " (+ mingentime, optimumT)"),
          ("pt_recognise", "reCOGnise specI cluster"), ("pt_", "porTraits: tool|trait[|model]; _p probability, _b binary"),
          ("rgi_", "RGI (CARD)"), ("bacmet_", "BacMet2 EXP"), ("as_", "antiSMASH"), ("gc_", "GECCO"), ("bs_", "BiG-SLiCE GCFs")]
with open(out.replace(".tsv.gz", ".columns.tsv"), "w") as fh:
    fh.write("column\tgroup\n")
    for c in M.columns:
        fh.write("{}\t{}\n".format(c, next((g for pre, g in groups if c.startswith(pre)), "")))
