#!/usr/bin/env python3
"""
Check and collate per-sample CoverM genome-mode outputs (one *_output_coverm.tsv per sample)
against the dereplicated MAG set and the Skyline sample sheet.

Outputs (results/coverm/):
  coverm_check_report.txt          everything printed below
  coverm_long.tsv.gz               sample x genome, all metrics, tidy format
  coverm_<metric>_matrix.tsv       genome x sample wide matrices
  coverm_rel_abund_masked_cf<T>.tsv  rel. abundance set to 0 where covered fraction < T
  coverm_sample_summary.tsv        per-sample QC + matched metadata
  coverm_genome_summary.tsv        per-genome prevalence / abundance
"""
import os, re, glob, argparse
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--coverm-dir", default=os.environ.get("COVERM_EXISTING"))
p.add_argument("--files", nargs="+", help="explicit list of CoverM tables (overrides --coverm-dir)")
p.add_argument("--metadata",   default=os.environ.get("METADATA"))
p.add_argument("--genomes",    default=os.environ.get("GENOMES"))
p.add_argument("--ext",        default=os.environ.get("GENOME_EXT", "fa"))
p.add_argument("--out",        default=os.path.join(os.environ.get("RESULTS", "."), "coverm"))
p.add_argument("--min-covfrac",  type=float, default=float(os.environ.get("COVERM_MIN_COVFRAC", 0.1)))
p.add_argument("--max-unmapped", type=float, default=50.0, help="flag samples above this unmapped %%")
a = p.parse_args()
os.makedirs(a.out, exist_ok=True)

REPORT = []
def log(s=""):
    print(s); REPORT.append(str(s))

METRICS = {"Mean": "mean", "Trimmed Mean": "trimmed_mean", "Relative Abundance (%)": "rel_abund",
           "Covered Fraction": "covered_fraction", "Read Count": "read_count",
           "RPKM": "rpkm", "TPM": "tpm", "Covered Bases": "covered_bases", "Length": "length"}

# ---------------- 1. per-sample CoverM files ----------------
files = sorted(a.files) if a.files else sorted(glob.glob(os.path.join(a.coverm_dir, "*_output_coverm.tsv")))
expected = {os.path.basename(f)[: -(len(a.ext) + 1)] for f in glob.glob(os.path.join(a.genomes, f"*.{a.ext}"))}
log(f"CoverM files: {len(files)} | dereplicated genomes: {len(expected)}")

frames, unmapped, bad_label, genome_issues, multi = [], {}, [], [], []
for f in files:
    sid = os.path.basename(f).replace("_output_coverm.tsv", "")
    df = pd.read_csv(f, sep="\t", na_values=["NA"])
    gcol, ren, labels = df.columns[0], {}, set()
    # longest metric names first so "Trimmed Mean" isn't caught by "Mean"
    for c in df.columns[1:]:
        for m in sorted(METRICS, key=len, reverse=True):
            if c.endswith(" " + m):
                ren[c] = METRICS[m]; labels.add(c[: -len(m) - 1]); break
    if len(labels) != 1:
        multi.append((sid, sorted(labels))); continue
    label = labels.pop()
    if not re.match(rf"^{re.escape(sid)}([_.]|$)", label):
        bad_label.append((sid, label))
    df = df.rename(columns={**ren, gcol: "genome"})
    um = df[df.genome == "unmapped"]
    unmapped[sid] = um["rel_abund"].iloc[0] if len(um) and "rel_abund" in um else np.nan
    df = df[df.genome != "unmapped"]
    found = set(df.genome)
    if found != expected or df.genome.duplicated().any():
        genome_issues.append((sid, len(found), len(expected - found), len(found - expected)))
    df.insert(0, "sample", sid)
    df["read_label"] = label
    frames.append(df)

if not frames:
    raise SystemExit(f"No readable *_output_coverm.tsv files in {a.coverm_dir}")
long = pd.concat(frames, ignore_index=True)
metrics = [m for m in METRICS.values() if m in long.columns]
log(f"Metrics present: {metrics}")

log("\n== File integrity ==")
log(f"Files with more/fewer than one sample in them: {len(multi)} {multi[:5]}")
log(f"Files whose column label doesn't match the file name: {len(bad_label)} {bad_label[:5]}")
log(f"Files whose genome set != dereplicated set: {len(genome_issues)}")
for s in genome_issues[:10]:
    log(f"   {s[0]}: {s[1]} genomes, {s[2]} missing, {s[3]} extra")

# ---------------- 2. per-sample QC ----------------
thr = a.min_covfrac
g = long.groupby("sample")
ss = pd.DataFrame({
    "read_label": g.read_label.first(),
    "unmapped_pct": pd.Series(unmapped),
    "genomes_rel_abund_sum": g.rel_abund.sum() if "rel_abund" in long else np.nan,
    "mapped_reads": g.read_count.sum() if "read_count" in long else np.nan,
    f"n_genomes_cf>={thr}": g.covered_fraction.apply(lambda x: (x >= thr).sum()) if "covered_fraction" in long else np.nan,
})
ss["total_pct"] = ss.genomes_rel_abund_sum + ss.unmapped_pct

log("\n== Per-sample ==")
log(f"Unmapped %: median {ss.unmapped_pct.median():.1f}, IQR {ss.unmapped_pct.quantile(.25):.1f}-"
    f"{ss.unmapped_pct.quantile(.75):.1f}, range {ss.unmapped_pct.min():.1f}-{ss.unmapped_pct.max():.1f}")
hi = ss[ss.unmapped_pct > a.max_unmapped].sort_values("unmapped_pct", ascending=False)
log(f"Samples > {a.max_unmapped}% unmapped: {len(hi)}")
if len(hi): log(hi[["unmapped_pct"]].head(15).to_string())
off = ss[(ss.total_pct - 100).abs() > 1]
log(f"Samples where genomes + unmapped != 100% (±1): {len(off)} {list(off.index[:5])}")
col = f"n_genomes_cf>={thr}"
log(f"Genomes detected per sample (covered fraction >= {thr}): median {ss[col].median():.0f}, "
    f"range {ss[col].min():.0f}-{ss[col].max():.0f}")
lo = ss.sort_values("mapped_reads").head(5)
log("Lowest mapped-read samples:\n" + lo[["mapped_reads", "unmapped_pct"]].to_string())

# ---------------- 3. per-genome QC ----------------
gg = long.groupby("genome")
gs = pd.DataFrame({
    "prevalence_cf": gg.covered_fraction.apply(lambda x: (x >= thr).sum()),
    "max_covered_fraction": gg.covered_fraction.max(),
    "mean_rel_abund": gg.rel_abund.mean(),
    "max_rel_abund": gg.rel_abund.max(),
    "total_reads": gg.read_count.sum() if "read_count" in long else np.nan,
}).sort_values("mean_rel_abund", ascending=False)
n_samp = long["sample"].nunique()
log("\n== Per-genome ==")
log(f"Genomes never reaching covered fraction {thr} in any sample: {(gs.prevalence_cf == 0).sum()}")
log(f"Genomes with zero reads in every sample: {(gs.total_reads == 0).sum()}")
log(f"Genomes detected in >= 50% of samples: {(gs.prevalence_cf >= n_samp / 2).sum()}")
log("Top 10 by mean relative abundance:\n" + gs.head(10)[["prevalence_cf", "mean_rel_abund"]].to_string())

# ---------------- 4. metadata ----------------
def read_meta(path):
    with open(path) as fh:
        header = [h.strip().replace(" ", "") for h in fh.readline().rstrip("\n").split("\t")]
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    width = max(len(r) for r in rows)
    if width == len(header) + 1 and "dna_plate_numdna_plate_col" in header:
        i = header.index("dna_plate_numdna_plate_col")
        header[i:i + 1] = ["dna_plate_num", "dna_plate_col"]
        log("NOTE: metadata header had 'dna_plate_numdna_plate_col' fused (missing tab) -- split it.")
    if width != len(header):
        log(f"WARNING: metadata rows have {width} fields but header has {len(header)} -- check column alignment!")
    rows = [(r + [""] * len(header))[: len(header)] for r in rows]
    return pd.DataFrame(rows, columns=header).apply(lambda s: s.str.strip())

log("\n== Metadata ==")
meta = read_meta(a.metadata)
by_tube = {v: i for i, v in meta.dna_tube_num.items() if v}
by_code = {v: i for i, v in meta.ceh_sample_code.items() if v}
match, mtype = {}, {}
for s in ss.index:
    if s in by_tube:   match[s], mtype[s] = by_tube[s], "dna_tube_num"
    elif s in by_code: match[s], mtype[s] = by_code[s], "ceh_sample_code"
unmatched = sorted(set(ss.index) - set(match))
log(f"Matched by dna_tube_num: {sum(v == 'dna_tube_num' for v in mtype.values())} | "
    f"by ceh_sample_code: {sum(v == 'ceh_sample_code' for v in mtype.values())} | unmatched: {len(unmatched)}")
if unmatched: log(f"   unmatched samples: {unmatched}")

dup = pd.Series(match).loc[lambda x: x.duplicated(keep=False)]
if len(dup):
    log("Same metadata row matched by >1 CoverM file (same biological sample sequenced twice?):")
    for row, grp in dup.groupby(dup):
        log(f"   {meta.loc[row, 'ceh_sample_code']}: {list(grp.index)}")
no_cov = meta[(meta.dna_tube_num != "") & ~meta.index.isin(match.values())]
log(f"Metadata rows with a DNA tube but no CoverM file: {len(no_cov)} {list(no_cov.dna_tube_num[:10])}")

m = meta.loc[[match[s] for s in ss.index if s in match]].copy()
m.index = [s for s in ss.index if s in match]
ss = ss.join(m.add_prefix("meta_")).assign(meta_match=pd.Series(mtype))

# ---------------- 5. write ----------------
long.to_csv(os.path.join(a.out, "coverm_long.tsv.gz"), sep="\t", index=False)
for met in metrics:
    long.pivot(index="genome", columns="sample", values=met).to_csv(
        os.path.join(a.out, f"coverm_{met}_matrix.tsv"), sep="\t")
if {"rel_abund", "covered_fraction"} <= set(long.columns):
    masked = long.assign(rel_abund=np.where(long.covered_fraction >= thr, long.rel_abund, 0))
    masked.pivot(index="genome", columns="sample", values="rel_abund").to_csv(
        os.path.join(a.out, f"coverm_rel_abund_masked_cf{thr}.tsv"), sep="\t")
ss.to_csv(os.path.join(a.out, "coverm_sample_summary.tsv"), sep="\t", index_label="sample")
gs.to_csv(os.path.join(a.out, "coverm_genome_summary.tsv"), sep="\t", index_label="genome")
with open(os.path.join(a.out, "coverm_check_report.txt"), "w") as fh:
    fh.write("\n".join(REPORT) + "\n")
log(f"\nWritten to {a.out}")
