#!/usr/bin/env python3
"""
Build the contig -> MAG -> dMAG table used by dereplicated_mag_link from a dRep collection.

  Cdb.csv  : every genome and its secondary cluster
  Wdb.csv  : the winner (dereplicated MAG, dMAG) of each cluster
  MAG FASTA: contig names, read from --genome-dirs (genomes whose FASTA is not found are reported, not linked)

Each MAG is assigned to the assembly it was binned from via its name prefix (text before "_bin."):
  --prefix-map Bulk_hifi_co=hifi_coassembly_bulk ...   explicit prefixes
  --default-assembly 'hifi_individual_{sample}'        otherwise, with {sample} = prefix minus Bulk_/Rhizosphere_

Safety check: for every prefix, the first contig of one MAG must be found, with an identical sequence, in that
assembly's contigs file (--assembly-contigs '<dir>/{asm}/contigs/contigs.fa'). metaMDBG numbers contigs per run,
so a MAG binned from a different assembly version would otherwise be linked silently to the wrong contigs.

Output (TSV): contig, asm, MAG, dMAG   (genome names without the .fa extension, as in GTDB-Tk user_genome)
"""
import argparse
import csv
import os
import sys
import time


def log(msg):
    print(time.strftime("[%Y-%m-%d %H:%M:%S] ") + msg, flush=True)


def strip_ext(name):
    for ext in (".fa", ".fna", ".fasta"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def prefix_of(genome):
    return genome.split("_bin.")[0] if "_bin." in genome else strip_ext(genome)


def read_fasta_records(path, wanted=None):
    """yield (name, sequence) - only for names in `wanted` if given"""
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None and (wanted is None or name in wanted):
        yield name, "".join(seq)


def fetch_contig(path, contig):
    """sequence of one contig from a (possibly large) FASTA; uses the .fai index if present"""
    fai = path + ".fai"
    if os.path.exists(fai):
        with open(fai) as f:
            for line in f:
                c, length, offset, linebases, linewidth = line.rstrip("\n").split("\t")[:5]
                if c == contig:
                    length, offset, linebases, linewidth = int(length), int(offset), int(linebases), int(linewidth)
                    nlines = (length + linebases - 1) // linebases
                    with open(path) as fa:
                        fa.seek(offset)
                        raw = fa.read(nlines * linewidth)
                    return raw.replace("\n", "").replace("\r", "")[:length]
        return None
    for _, s in read_fasta_records(path, wanted={contig}):
        return s
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cdb", required=True)
    ap.add_argument("--wdb", required=True)
    ap.add_argument("--genome-dirs", nargs="+", required=True)
    ap.add_argument("--prefix-map", nargs="*", default=[])
    ap.add_argument("--default-assembly", default="hifi_individual_{sample}")
    ap.add_argument("--assembly-contigs", required=True, help="path template with {asm}")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    pmap = dict(x.split("=", 1) for x in a.prefix_map if x)

    def assembly_of(genome):
        p = prefix_of(genome)
        if p in pmap:
            return pmap[p]
        sample = p
        for comp in ("Bulk_", "Rhizosphere_"):
            if sample.startswith(comp):
                sample = sample[len(comp):]
        return a.default_assembly.format(sample=sample)

    with open(a.wdb) as f:
        winner = {r["cluster"]: r["genome"] for r in csv.DictReader(f)}
    with open(a.cdb) as f:
        cdb = [(r["genome"], r["secondary_cluster"]) for r in csv.DictReader(f)]
    log("{} genomes in {} clusters ({} winners)".format(len(cdb), len({c for _, c in cdb}), len(winner)))

    # locate MAG FASTA files
    found, missing = {}, []
    for genome, _ in cdb:
        for d in a.genome_dirs:
            path = os.path.join(d, genome)
            if os.path.isfile(path) and os.access(path, os.R_OK):   # exists AND readable (permissions differ per folder)
                found[genome] = path
                break
        else:
            missing.append(genome)
    unreadable = [g for g in missing if any(os.path.exists(os.path.join(d, g)) for d in a.genome_dirs)]
    log("{} MAG FASTA files readable, {} missing or unreadable ({} exist but lack read permission); "
        "their contigs cannot be linked".format(len(found), len(missing), len(unreadable)))
    winners_missing = [g for g in winner.values() if g not in found]
    if winners_missing:
        sys.exit("dereplicated MAGs (winners) missing from --genome-dirs, e.g. {}".format(winners_missing[:3]))

    # safety check: one contig per prefix must match its assembly exactly
    checked = {}
    for genome, path in sorted(found.items()):
        p = prefix_of(genome)
        if p in checked:
            continue
        asm = assembly_of(genome)
        asm_path = a.assembly_contigs.format(asm=asm)
        if not os.path.exists(asm_path):
            sys.exit("assembly contigs not found for prefix {} ({}): {}".format(p, asm, asm_path))
        contig, seq = next(read_fasta_records(path))
        ref = fetch_contig(asm_path, contig)
        if ref is None or ref.upper() != seq.upper():
            sys.exit("prefix {}: contig {} of {} is {} in assembly {} ({}). The MAGs were not binned from this "
                     "assembly - fix dmag_link.mag_prefix_to_assembly / default_assembly."
                     .format(p, contig, genome, "missing" if ref is None else "different", asm, asm_path))
        checked[p] = asm
        log("  ok: {:<30} -> {} (contig {}, {} bp identical)".format(p, asm, contig, len(seq)))

    # write the table
    tmp = a.out + ".part"
    n = 0
    with open(tmp, "w") as out:
        out.write("contig\tasm\tMAG\tdMAG\n")
        for genome, cluster in cdb:
            if genome not in found:
                continue
            asm = assembly_of(genome)
            mag, dmag = strip_ext(genome), strip_ext(winner[cluster])
            with open(found[genome]) as f:
                for line in f:
                    if line.startswith(">"):
                        out.write("{}\t{}\t{}\t{}\n".format(line[1:].split()[0], asm, mag, dmag))
                        n += 1
    os.replace(tmp, a.out)
    log("wrote {} contigs from {} MAGs ({} assemblies) to {}".format(n, len(found), len(checked), a.out))


if __name__ == "__main__":
    main()
