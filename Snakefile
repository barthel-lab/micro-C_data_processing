include: "workflow/pipeline.smk"

# Define QC output files of workflows
rule all:
   input:
        expand("results/finalBam/{aliquot_barcode}.mapped.PT.bam", aliquot_barcode=Sname),
        expand("results/finalBam/{aliquot_barcode}.mapped.PT.bam.bai", aliquot_barcode=Sname),
        expand("results/LibQC/{aliquot_barcode}.QC.txt", aliquot_barcode=Sname),
        expand("results/libComX/{aliquot_barcode}.libComX.preseq", aliquot_barcode=Sname),
        expand("results/contacMatrix/{aliquot_barcode}.contact_map.hic", aliquot_barcode=Sname),
        expand("results/coolerMatrix/{aliquot_barcode}.matrix_1kb.cool", aliquot_barcode=Sname),
        expand("results/coolerMatrix/{aliquot_barcode}.matrix_1kb.mcool", aliquot_barcode=Sname)

