
# General micro-C pipeline
#include: "workflow/pipeline.smk"
# Special version to Telomere-C project
include: "workflow/pipeline2.smk"

# Define QC output files of workflows
rule all:
   input:
# General micro-C
      #  expand("results/LibQC/{aliquot_barcode}.QC.txt", aliquot_barcode=Sname),
      #  expand("results/bamCoverage/{aliquot_barcode}.filt.mapped.PT.bigwig", aliquot_barcode=Sname),
      #  expand("results/bamCoverage/{aliquot_barcode}.filt.unmapped.PT.bigwig", aliquot_barcode=Sname),
      #  expand("results/contacMatrix/{aliquot_barcode}.contact_map.hic", aliquot_barcode=Sname),
      #  expand("results/coolerMatrix/{aliquot_barcode}.matrix_1kb.cool", aliquot_barcode=Sname),
      #  expand("results/coolerMatrix/{aliquot_barcode}.matrix_1kb.mcool", aliquot_barcode=Sname),
      #  expand("results/plotCircos/{aliquot_barcode}.circos.png", aliquot_barcode=Sname)
# Telomere-C special
        expand("results/LibQC_filt/{aliquot_barcode}.filt.QC.txt", aliquot_barcode=Sname),
        expand("results/bamCoverage_filt/{aliquot_barcode}.filt.mapped.PT.bigwig", aliquot_barcode=Sname),
        expand("results/bamCoverage_filt/{aliquot_barcode}.filt.unmapped.PT.bigwig", aliquot_barcode=Sname),
        expand("results/contacMatrix_filt/{aliquot_barcode}.filt.contact_map.hic", aliquot_barcode=Sname),
        expand("results/coolerMatrix_filt/{aliquot_barcode}.filt.matrix_1kb.cool", aliquot_barcode=Sname),
        expand("results/coolerMatrix_filt/{aliquot_barcode}.filt.matrix_1kb.mcool", aliquot_barcode=Sname),
        expand("results/plotCircos_filt/{aliquot_barcode}.circos.png", aliquot_barcode=Sname)
        

