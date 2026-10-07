#!/usr/bin/env python3
"""Make GTDB-Tk optional in porTraits' main.nf (params.skip_gtdbtk). Idempotent."""
import sys
p = sys.argv[1]
s = open(p).read()
if "params.skip_gtdbtk" in s:
    sys.exit(0)
old = ("\t// GTDBtk classification\t\n"
       "\tgtdbtk_classify(genomes_ch, params.gtdbtk_data)\n\n"
       "\tall_results_ch = all_results_ch\n"
       "\t\t.mix(gtdbtk_classify.out.gtdb_taxonomy)\n")
new = ("\t// GTDBtk classification (optional; skyline patch)\n"
       "\tif (!params.skip_gtdbtk) {\n"
       "\t\tgtdbtk_classify(genomes_ch, params.gtdbtk_data)\n\n"
       "\t\tall_results_ch = all_results_ch\n"
       "\t\t\t.mix(gtdbtk_classify.out.gtdb_taxonomy)\n"
       "\t}\n")
if old not in s:
    sys.exit("GTDB-Tk block not found in " + p + " -- porTraits version changed?")
open(p, "w").write(s.replace(old, new))
print("patched", p)
