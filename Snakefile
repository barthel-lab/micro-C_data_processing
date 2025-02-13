include: "workflow/pipeline.smk"
include: "workflow/pipeline2.smk"

# Define QC output files of workflows
rule all:
   input:
        expand("results/finalBam/{aliquot_barcode}.mapped.PT.bam", aliquot_barcode=Sname),
        expand("results/finalBam/{aliquot_barcode}.mapped.PT.bam.bai", aliquot_barcode=Sname),
        expand("results/LibQC/{aliquot_barcode}.QC.txt", aliquot_barcode=Sname),
        expand("results/libComX/{aliquot_barcode}.libComX.preseq", aliquot_barcode=Sname),
        expand("results/contacMatrix/{aliquot_barcode}.contact_map.hic", aliquot_barcode=Sname),
        expand("results/coolerMatrix/{aliquot_barcode}.matrix_1kb.cool", aliquot_barcode=Sname),
        expand("results/coolerMatrix/{aliquot_barcode}.matrix_1kb.mcool", aliquot_barcode=Sname),
        expand("results/KMCdb/{aliquot_barcode}.{read}.{sub}", aliquot_barcode=Sname, read=["R1","R2"], sub=["kmc_pre","kmc_suf"]),
        expand("results/KMCDump/{aliquot_barcode}.{read}.dump.txt", aliquot_barcode=Sname, read=["R1","R2"]),
        expand("results/KMCDump/{aliquot_barcode}.{read}.dump.subset.txt", aliquot_barcode=Sname, read=["R1","R2"]),
        expand("results/patternList/{aliquot_barcode}.{read}.dump.pattern.txt", aliquot_barcode=Sname, read=["R1","R2"]),
        expand("results/preFastq/{aliquot_barcode}.{read}.filt.fastq.gz", aliquot_barcode=Sname, read=["R1","R2"]),
        expand("results/preFastq/{aliquot_barcode}.readname.txt", aliquot_barcode=Sname),
        expand("results/preFastq/{aliquot_barcode}.readname.tmp.txt", aliquot_barcode=Sname),
        expand("results/finalBam_filt/{aliquot_barcode}.filt.mapped.PT.bam", aliquot_barcode=Sname),
        expand("results/finalBam_filt/{aliquot_barcode}.filt.mapped.PT.bam.bai", aliquot_barcode=Sname),
        expand("results/LibQC_filt/{aliquot_barcode}.filt.QC.txt", aliquot_barcode=Sname),
        expand("results/libComX_filt/{aliquot_barcode}.filt.libComX.preseq", aliquot_barcode=Sname),
        expand("results/contacMatrix_filt/{aliquot_barcode}.filt.contact_map.hic", aliquot_barcode=Sname),
        expand("results/coolerMatrix_filt/{aliquot_barcode}.filt.matrix_1kb.cool", aliquot_barcode=Sname),
        expand("results/coolerMatrix_filt/{aliquot_barcode}.filt.matrix_1kb.mcool", aliquot_barcode=Sname)
 

