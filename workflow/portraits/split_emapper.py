#!/usr/bin/env python3
"""Split one eggNOG-mapper .emapper.annotations file (queries named "<mag>|<protein>") into
<outdir>/<mag>/<mag>.emapper.annotations, one per MAG listed in --mags (header-only file if a MAG has no hits),
with the "<mag>|" prefix removed so each file looks like a per-genome eggNOG-mapper run."""
import argparse, os, sys

p = argparse.ArgumentParser()
p.add_argument("--annotations", required=True)
p.add_argument("--mags", nargs="+", required=True)
p.add_argument("--outdir", required=True)
a = p.parse_args()

prelude, header, rows, trailer = [], None, {m: [] for m in a.mags}, []
with open(a.annotations) as fh:
    for line in fh:
        if line.startswith("##"):
            (prelude if header is None else trailer).append(line)
        elif line.startswith("#query"):
            header = line
        else:
            mag, sep, rest = line.partition("|")
            if not sep or mag not in rows:
                sys.exit("unexpected query id (no known '<mag>|' prefix): " + line.split("\t")[0])
            rows[mag].append(rest)
if header is None:
    sys.exit("no '#query' header found in " + a.annotations)

for mag, lines in rows.items():
    d = os.path.join(a.outdir, mag); os.makedirs(d, exist_ok=True)
    tmp = os.path.join(d, mag + ".emapper.annotations.tmp")
    with open(tmp, "w") as out:
        out.writelines(prelude + [header] + lines + trailer)
    os.replace(tmp, os.path.join(d, mag + ".emapper.annotations"))
    print(mag, len(lines))
