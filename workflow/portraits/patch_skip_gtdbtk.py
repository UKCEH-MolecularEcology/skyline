#!/usr/bin/env python3
"""skyline patches for porTraits (idempotent). Usage: patch_skip_gtdbtk.py <porTraits>/main.nf
 1. main.nf: GTDB-Tk optional (params.skip_gtdbtk)
 2. modules/recognise.nf: MAGs without reCOGnise marker genes get an empty cogs.txt
    (recognise exits 0 but omits it; proteins/genes/gff are written, so all trait predictors still run)
 3. main.nf: optional precomputed eggNOG-mapper annotations (params.emapper_annotations = directory with
    <id>/<id>.emapper.annotations or <id>.emapper.annotations); eggNOG-mapper only runs for genomes without one
"""
import os, sys

main_nf = sys.argv[1]
rec_nf = os.path.join(os.path.dirname(main_nf), "portraits", "modules", "recognise.nf")

def patch(path, old, new, tag):
    s = open(path).read()
    if tag in s:
        print("already patched:", path); return
    if old not in s:
        sys.exit("pattern not found in {} -- porTraits version changed?".format(path))
    open(path, "w").write(s.replace(old, new)); print("patched:", path)

patch(main_nf,
      ("\t// GTDBtk classification\t\n"
       "\tgtdbtk_classify(genomes_ch, params.gtdbtk_data)\n\n"
       "\tall_results_ch = all_results_ch\n"
       "\t\t.mix(gtdbtk_classify.out.gtdb_taxonomy)\n"),
      ("\t// GTDBtk classification (optional; skyline patch)\n"
       "\tif (!params.skip_gtdbtk) {\n"
       "\t\tgtdbtk_classify(genomes_ch, params.gtdbtk_data)\n\n"
       "\t\tall_results_ch = all_results_ch\n"
       "\t\t\t.mix(gtdbtk_classify.out.gtdb_taxonomy)\n"
       "\t}\n"),
      "params.skip_gtdbtk")

patch(rec_nf,
      "\trm -fv genome_file\n",
      ("\t# skyline patch: no marker genes found -> recognise writes no cogs.txt\n"
       "\t[ -e ${genome_id}/recognise/${genome_id}.cogs.txt ] || touch ${genome_id}/recognise/${genome_id}.cogs.txt\n"
       "\trm -fv genome_file\n"),
      "skyline patch: no marker genes")

patch(main_nf,
      ("\t// emapper\n"
       "\teggnog_mapper(\n"
       "\t\trecognise_genome.out.proteins,\n"
       "\t\tparams.eggnog_db\n"
       "\t)\n\n"
       "\temapper2matrix(\n"
       "\t\teggnog_mapper.out.eggnog,\n"
       "\t\tparams.pfam_clade_map\n"
       "\t)\n"),
      ("\t// emapper (skyline patch: reuse precomputed annotations where available)\n"
       "\tdef emapper_file = { gid ->\n"
       "\t\tif (!params.emapper_annotations) return null\n"
       "\t\tdef nested = file(\"${params.emapper_annotations}/${gid}/${gid}.emapper.annotations\")\n"
       "\t\tdef flat = file(\"${params.emapper_annotations}/${gid}.emapper.annotations\")\n"
       "\t\treturn nested.exists() ? nested : (flat.exists() ? flat : null)\n"
       "\t}\n"
       "\temapper_in_ch = recognise_genome.out.proteins\n"
       "\t\t.branch { gid, proteins ->\n"
       "\t\t\tprecomputed: emapper_file(gid) != null\n"
       "\t\t\tmissing: true\n"
       "\t\t}\n\n"
       "\teggnog_mapper(\n"
       "\t\temapper_in_ch.missing,\n"
       "\t\tparams.eggnog_db\n"
       "\t)\n\n"
       "\teggnog_ch = emapper_in_ch.precomputed\n"
       "\t\t.map { gid, proteins -> tuple(gid, emapper_file(gid)) }\n"
       "\t\t.mix(eggnog_mapper.out.eggnog)\n\n"
       "\temapper2matrix(\n"
       "\t\teggnog_ch,\n"
       "\t\tparams.pfam_clade_map\n"
       "\t)\n"),
      "skyline patch: reuse precomputed annotations")
