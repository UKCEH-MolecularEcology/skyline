# scripts/analyse_uc_clusters.py

import pandas as pd
from collections import defaultdict

uc_file = snakemake.input.uc
out_mixed = snakemake.output.mixed
out_amp = snakemake.output.amp_only
out_lr = snakemake.output.lr_only
out_sing_amp = snakemake.output.sing_amp
out_sing_lr = snakemake.output.sing_lr
log_file = snakemake.log[0]

def classify_source(seq_id):
    return "amplicon" if seq_id.startswith("ASV") else "longread"

def parse_uc_file(filepath):
    cluster_members = defaultdict(list)
    cluster_reps = {}

    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            record_type = parts[0]
            cluster_num = parts[1]
            seq_id = parts[8].replace("^M", "")
            if record_type == "S":
                cluster_reps[cluster_num] = seq_id
                cluster_members[cluster_num].append(seq_id)
            elif record_type == "H":
                cluster_members[cluster_num].append(seq_id)
    return cluster_members, cluster_reps

def summarize_clusters(cluster_members, cluster_reps):
    amplicon_only = []
    longread_only = []
    mixed = []
    singleton_amplicon = []
    singleton_longread = []

    for cluster, members in cluster_members.items():
        sources = [classify_source(m) for m in members]
        rep = cluster_reps.get(cluster, "NA")
        other_members = [m for m in members if m != rep]
        members_str = ";".join(other_members)

        if len(members) == 1:
            if sources[0] == "amplicon":
                singleton_amplicon.append(members[0])
            else:
                singleton_longread.append(members[0])
        elif all(s == "amplicon" for s in sources):
            amplicon_only.append((cluster, rep, len(members), members_str))
        elif all(s == "longread" for s in sources):
            longread_only.append((cluster, rep, len(members), members_str))
        else:
            amplicon_count = sources.count("amplicon")
            longread_count = sources.count("longread")
            mixed.append((cluster, rep, amplicon_count, longread_count, members_str))

    return amplicon_only, longread_only, mixed, singleton_amplicon, singleton_longread

clusters, reps = parse_uc_file(uc_file)
amp_only, lr_only, mixed, sing_amp, sing_lr = summarize_clusters(clusters, reps)

pd.DataFrame(mixed, columns=["Cluster", "Representative", "Amplicon_Count", "Longread_Count", "Members"]).to_csv(out_mixed, index=False)
pd.DataFrame(amp_only, columns=["Cluster", "Representative", "Size", "Members"]).to_csv(out_amp, index=False)
pd.DataFrame(lr_only, columns=["Cluster", "Representative", "Size", "Members"]).to_csv(out_lr, index=False)
pd.DataFrame(sing_amp, columns=["Sequence"]).to_csv(out_sing_amp, index=False)
pd.DataFrame(sing_lr, columns=["Sequence"]).to_csv(out_sing_lr, index=False)

with open(log_file, "w") as log:
    log.write(f"Amplicon-only clusters: {len(amp_only)}\n")
    log.write(f"Long-read-only clusters: {len(lr_only)}\n")
    log.write(f"Mixed clusters: {len(mixed)}\n")
    log.write(f"Amplicon singletons: {len(sing_amp)}\n")
    log.write(f"Long-read singletons: {len(sing_lr)}\n")
