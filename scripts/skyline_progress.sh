#!/bin/bash
# Progress of the MAG annotation run (porTraits, antiSMASH, GECCO, BiG-SLiCE).
# Paths, MAG count and batch size are read from config/config.yaml.
# Usage (from the repo root, on the head node):
#   bash scripts/skyline_progress.sh
#   watch -n 300 bash scripts/skyline_progress.sh

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CFG="${REPO}/config/config.yaml"
cfg() { grep -m1 "^  $1:" "$CFG" | sed "s/^  $1: *//; s/\"//g; s/ *#.*//"; }   # keys under annotations:

A=$(cfg results_dir)
MAGS_DIR=$(cfg mags_dir); EXT=$(cfg mags_ext); EXT=${EXT:-fa}
BS=$(cfg batch_size); BS=${BS:-100}
N_MAG=$(ls "$MAGS_DIR"/*."$EXT" 2>/dev/null | wc -l)
N_BATCH=$(( (N_MAG + BS - 1) / BS ))
count() { ls $1 2>/dev/null | wc -l; }

echo "== $(date '+%a %d %b %H:%M')  ($N_MAG MAGs, $N_BATCH batches of $BS)"
echo "-- batches done"
printf "  %-10s %4s / %s\n" porTraits "$(count "$A/portraits/batches/*/collated/portraits_results.tsv.gz")" $N_BATCH
printf "  %-10s %4s / %s\n" antiSMASH "$(count "$A/bgc/antismash/batches/*.done")" $N_BATCH
printf "  %-10s %4s / %s\n" GECCO     "$(count "$A/bgc/gecco/batches/*.done")" $N_BATCH

echo "-- MAGs"
ok=$(grep -l '^ok' $A/bgc/antismash/per_mag/*/STATUS 2>/dev/null | wc -l)
fail=$(grep -l '^failed' $A/bgc/antismash/per_mag/*/STATUS 2>/dev/null | wc -l)
printf "  %-10s %5s / %s ok, %s failed\n" antiSMASH $ok $N_MAG $fail
gdone=$(ls -d $A/bgc/gecco/per_mag/*/ 2>/dev/null | wc -l)
gfail=$(count "$A/bgc/gecco/per_mag/*.FAILED")
printf "  %-10s %5s / %s done, %s failed\n" GECCO $gdone $N_MAG $gfail

echo "-- porTraits skipped tasks (finished batches)"
grep -h "skipped tasks" $A/logs/portraits/batch*.log 2>/dev/null | sed 's/(.*//' | sort | uniq -c | sed 's/^/  /'

echo "-- final outputs"
for f in portraits/portraits_results.tsv.gz bgc/antismash_regions_all.tsv bgc/gecco_clusters_all.tsv \
         bgc/bgc_overlap.tsv bgc/bigslice/tables/gcf_summary.tsv; do
  if [ -e $A/$f ]; then printf "  %-45s %s\n" $f "$(date -r $A/$f '+%d %b %H:%M')"
  else printf "  %-45s -\n" $f; fi
done

echo "-- SLURM jobs"
if command -v squeue >/dev/null; then squeue -u "$USER" -h -o "%j %T" | sort | uniq -c | sed 's/^/  /'
else echo "  (squeue not available here -- run on the head node)"; fi
