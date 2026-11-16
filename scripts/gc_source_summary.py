#!/usr/bin/env python3
"""
Gene catalogue: number of clusters per level and where their genes come from.

Gene names start with their input's prefix (sample = <prefix>_<name>, gene = <sample>|<orf>):
  sr_  short-read assemblies     lr_  long-read assemblies     mag_  per-MAG ORF calls
Each cluster gets the set of sources of all genes it contains; levels are chained (100 -> 99 -> 95 -> ...),
so a level-L cluster's sources are the union over the level-(L+1) clusters it contains.

Output (tab-separated text): one row per level with the number of clusters and the number / percentage that
contain genes only from short-read assemblies, only from long-read sources (long-read assemblies and/or MAG ORFs),
from both (shared), and only from MAG ORFs (genes called only at MAG level).
"""
import argparse
import os
import time

import duckdb


def log(msg):
    print(time.strftime("[%Y-%m-%d %H:%M:%S] ") + msg, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", required=True, help="clustering levels in order, e.g. 100,99,95,90,70,50")
    ap.add_argument("--clusters-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--memory", default="300GB")
    a = ap.parse_args()
    levels = [int(x) for x in a.levels.split(",")]
    os.makedirs(a.tmp, exist_ok=True)

    con = duckdb.connect(os.path.join(a.tmp, "sources.duckdb"))
    con.execute("SET threads = {}".format(a.threads))
    con.execute("SET memory_limit = '{}'".format(a.memory))
    con.execute("SET temp_directory = '{}'".format(a.tmp))
    con.execute("SET preserve_insertion_order = false")
    clu = lambda lv: ("read_csv('{}', delim='\t', header=false, quote='', escape='', "
                      "columns={{'rep': 'VARCHAR', 'member': 'VARCHAR'}})").format(
                          os.path.join(a.clusters_dir, "clusters_{}.tsv.gz".format(lv)))
    # source bit per gene: 1 = short-read assembly, 2 = long-read assembly, 4 = MAG ORF
    bit = ("CASE WHEN starts_with(member, 'sr_') THEN 1 WHEN starts_with(member, 'lr_') THEN 2 "
           "WHEN starts_with(member, 'mag_') THEN 4 ELSE 8 END")

    rows = []
    prev = None
    for lv in levels:
        log("level {}".format(lv))
        if prev is None:   # first level: members are genes
            con.execute("CREATE OR REPLACE TABLE f_{} AS SELECT rep, bit_or({}) AS f FROM {} GROUP BY rep".format(lv, bit, clu(lv)))
        else:              # members are the previous level's representatives
            con.execute("CREATE OR REPLACE TABLE f_{} AS SELECT c.rep, bit_or(p.f) AS f FROM {} c "
                        "JOIN f_{} p ON c.member = p.rep GROUP BY c.rep".format(lv, clu(lv), prev))
            con.execute("DROP TABLE f_{}".format(prev))
        r = con.execute("""SELECT count(*),
                 sum(CASE WHEN f = 1 THEN 1 ELSE 0 END),
                 sum(CASE WHEN f & 1 = 0 AND f & 6 > 0 THEN 1 ELSE 0 END),
                 sum(CASE WHEN f & 1 = 1 AND f & 6 > 0 THEN 1 ELSE 0 END),
                 sum(CASE WHEN f = 4 THEN 1 ELSE 0 END),
                 sum(CASE WHEN f & 8 = 8 THEN 1 ELSE 0 END)
               FROM f_{}""".format(lv)).fetchone()
        rows.append((lv,) + tuple(int(x or 0) for x in r))
        log("  clusters {:,}".format(rows[-1][1]))
        prev = lv

    pct = lambda n, t: "{:.2f}".format(100.0 * n / t) if t else "NA"
    with open(a.out + ".part", "w") as out:
        out.write("# Gene catalogue clusters by source of their genes\n"
                  "# SR = short-read assemblies; LR = long-read assemblies and/or per-MAG ORF calls; "
                  "shared = both; MAG_only = genes found only in the per-MAG ORF calls (a subset of LR_only)\n")
        out.write("level\tclusters\tSR_only\tSR_only_pct\tLR_only\tLR_only_pct\tshared\tshared_pct\t"
                  "MAG_only\tMAG_only_pct\tunknown_source\n")
        for lv, n, sr, lr, sh, mag, unk in rows:
            out.write("\t".join(map(str, [lv, n, sr, pct(sr, n), lr, pct(lr, n), sh, pct(sh, n),
                                          mag, pct(mag, n), unk])) + "\n")
    os.replace(a.out + ".part", a.out)
    log("wrote " + a.out)


if __name__ == "__main__":
    main()
