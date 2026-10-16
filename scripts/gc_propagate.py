#!/usr/bin/env python3
"""
Propagate gene-catalogue annotations from the annotated level (e.g. 90) to another clustering level.

Levels are chained (100 -> 99 -> 95 -> 90 -> 70 -> 50): clusters_<L>.tsv.gz holds
(representative at L, member) where the members are the representatives of the level above.

  level < annotated (70, 50): every representative is itself an annotated representative
                              -> the annotation rows of those representatives
  level > annotated (95, 99, 100): each representative belongs to one annotated cluster
                              -> rep_<L> joined through the chain to its annotated representative,
                                 plus that representative's annotation (inherited, not re-annotated)

Outputs (gzip TSV, header): <tool>_reps.tsv.gz per tool; for levels above the annotated one also
reps_to_<annotated>.tsv.gz (rep_<L>, rep_<annotated>) for every representative, annotated or not.
"""
import argparse
import os
import sys
import time

import duckdb


def log(msg):
    print(time.strftime("[%Y-%m-%d %H:%M:%S] ") + msg, flush=True)


def read_clusters(path):
    return ("read_csv('{}', delim='\t', header=false, quote='', escape='', "
            "columns={{'rep': 'VARCHAR', 'member': 'VARCHAR'}})").format(path)


def read_table(path):
    return ("read_csv('{}', delim='\t', header=true, all_varchar=true, quote='', escape='', "
            "null_padding=true)").format(path)


def copy_out(con, query, path):
    tmp = path + ".part"
    con.execute("COPY ({}) TO '{}' (FORMAT csv, DELIMITER '\t', HEADER true, QUOTE '', "
                "COMPRESSION gzip)".format(query, tmp))
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--annot-level", type=int, required=True)
    ap.add_argument("--levels", required=True, help="clustering levels in order, e.g. 100,99,95,90,70,50")
    ap.add_argument("--clusters-dir", required=True)
    ap.add_argument("--annot-dir", required=True, help="annotation/<annot_level> directory")
    ap.add_argument("--tools", required=True, help="e.g. rgi,bacmet,eggnog")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--memory", default="100GB")
    a = ap.parse_args()

    levels = [int(x) for x in a.levels.split(",")]
    L, A = a.level, a.annot_level
    if L == A or L not in levels or A not in levels:
        sys.exit("level {} must differ from the annotated level {} and both be in {}".format(L, A, levels))
    os.makedirs(a.out_dir, exist_ok=True)
    os.makedirs(a.tmp, exist_ok=True)

    con = duckdb.connect()
    con.execute("SET threads = {}".format(a.threads))
    con.execute("SET memory_limit = '{}'".format(a.memory))
    con.execute("SET temp_directory = '{}'".format(a.tmp))
    con.execute("SET preserve_insertion_order = false")
    clu = lambda lv: os.path.join(a.clusters_dir, "clusters_{}.tsv.gz".format(lv))

    if L < A:
        # representatives at L are a subset of the annotated representatives
        log("level {}: distinct representatives from {}".format(L, clu(L)))
        con.execute("CREATE TEMP TABLE reps AS SELECT DISTINCT rep FROM {}".format(read_clusters(clu(L))))
        log("  {} representatives".format(con.execute("SELECT count(*) FROM reps").fetchone()[0]))
        keyed = "reps"
    else:
        # chain from L down to A: clusters_<l> for every level between them (exclusive of L, inclusive of A)
        chain = [lv for lv in levels if A <= lv < L]          # e.g. L=99, A=90 -> [95, 90]
        chain.sort(reverse=True)                               # 95, 90 (closest to L first)
        log("level {} -> {} via clusters {}".format(L, A, chain))
        # start: member of clusters_<first> is rep_<L>
        q = "SELECT member AS rep_l, rep AS cur FROM {}".format(read_clusters(clu(chain[0])))
        con.execute("CREATE TEMP TABLE m AS " + q)
        for lv in chain[1:]:
            log("  joining clusters_{}".format(lv))
            con.execute("CREATE TEMP TABLE m2 AS SELECT m.rep_l, c.rep AS cur FROM m "
                        "JOIN {} c ON c.member = m.cur".format(read_clusters(clu(lv))))
            con.execute("DROP TABLE m")
            con.execute("ALTER TABLE m2 RENAME TO m")
        n = con.execute("SELECT count(*), count(DISTINCT cur) FROM m").fetchone()
        log("  {} representatives at {} -> {} representatives at {}".format(n[0], L, n[1], A))
        out_map = os.path.join(a.out_dir, "reps_to_{}.tsv.gz".format(A))
        copy_out(con, "SELECT rep_l AS rep_{}, cur AS rep_{} FROM m".format(L, A), out_map)
        log("  wrote " + out_map)
        keyed = "m"

    for tool in a.tools.split(","):
        src = os.path.join(a.annot_dir, "{}_genes.tsv".format(tool))
        out = os.path.join(a.out_dir, "{}_reps.tsv.gz".format(tool))
        cols = [r[0] for r in con.execute("DESCRIBE SELECT * FROM {}".format(read_table(src))).fetchall()]
        key = cols[0]                                          # ORF_ID / gene / query = annotated representative
        log("{}: key column '{}' ({} columns)".format(tool, key, len(cols)))
        con.execute("CREATE OR REPLACE TEMP TABLE ann AS SELECT * FROM {}".format(read_table(src)))
        if L < A:
            q = 'SELECT a.* FROM ann a SEMI JOIN reps r ON a."{}" = r.rep'.format(key)
        else:
            q = ('SELECT m.rep_l AS rep_{}, a.* FROM m JOIN ann a ON a."{}" = m.cur'
                 ).format(L, key)
        copy_out(con, q, out)
        log("  wrote " + out)
        con.execute("DROP TABLE ann")

    log("done")


if __name__ == "__main__":
    main()
