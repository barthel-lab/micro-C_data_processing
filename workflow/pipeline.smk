# Configure
# Reference genome (Make sure coresponding index files are in the same director)
# Please refer this website of detailed alignment method https://micro-c.readthedocs.io/en/latest/fastq_to_bam.html
ref_fasta="/tgen_labs/barthel/references/CHM13v2/chm13v2.0.fasta"
ref_genome="/tgen_labs/barthel/references/CHM13v2/chm13v2.0.size.genome"
srcget_qc="/tgen_labs/barthel/software/github/external/Micro-C/get_qc.py"
srcJucier="/tgen_labs/barthel/software/github/external/Micro-C/juicer_tools_1.22.01.jar"

# Define sample names
import pandas as pd
fastqls = pd.read_csv("fastqList.txt", sep='\t', header=None, names=["name","R1","R2"])
fastqls.index = fastqls['name']
Sname = pd.Series(fastqls['name'])

# Align your Micro-C library to the reference. Please note the specific settings that are needed to map mates independently and for optimal results with our proximity library reads.
rule bwa:
    input:
        R1 = lambda wildcards: fastqls.loc[wildcards.aliquot_barcode].iloc[1],
        R2 = lambda wildcards: fastqls.loc[wildcards.aliquot_barcode].iloc[2]
    output:
        temp("results/bwa/{aliquot_barcode}.aln.sam")
    params:
        ref = ref_fasta
    threads: 32
    resources:
         mem_mb=128728,
         runtime=2159
    shell:"""
    bwa mem -5SP -T0 -t32 \
    {params.ref} \
    {input.R1} \
    {input.R2} \
    -o {output}"""

# We use the parse module of the pairtools pipeline to find ligation junctions in Micro-C (and other proximity ligation) libraries. When a ligation event is identified in the alignment file the pairtools pipeline will record the outer-most (5’) aligned base pair and the strand of each one of the paired reads into .pairsam file (pairsam format captures SAM entries together with the Hi-C pair information). In addition, it will also asign a pair type for each event. e.g. if both reads aligned uniquely to only one region in the genome, the type UU (Unique-Unique) will be assigned to the pair. The following steps are necessary to identify the high quality valid pairs over low quality events (e.g. due to low mapping quality):
rule RecordValidLigation:
    input:
        "results/bwa/{aliquot_barcode}.aln.sam"
    output:
        temp("results/RecordValidLigation/{aliquot_barcode}.parsed.pairsam")
    params:
        ref = ref_fasta
    conda:
        "micro-C"
    shell:"""
    pairtools parse \
    --min-mapq 40 \
    --walks-policy 5unique \
    --max-inter-align-gap 30 \
    --nproc-in 8 \
    --nproc-out 8 \
    --chroms-path {params.ref} \
    {input} > {output}"""

rule sortParisam:
    input:
        "results/RecordValidLigation/{aliquot_barcode}.parsed.pairsam"
    output:
        temp("results/sortParisam/{aliquot_barcode}.sort.pairsam")
    threads: 32
    resources:
        mem_mb=65536
    conda:
        "micro-C"
    shell:"""
    pairtools sort \
    --nproc 16 \
    --tmpdir=temp \
    {input} > {output}"""

rule rmPCRDup:
    input:
        "results/sortParisam/{aliquot_barcode}.sort.pairsam"
    output:
        pairsam = temp("results/rmPCRDup/{aliquot_barcode}.dedup.pairsam"),
        stats = "results/rmPCRDup/{aliquot_barcode}.stats.txt"
    threads: 16
    resources:
        mem_mb=32768
    conda:
        "micro-C"
    shell:"""
    pairtools dedup \
    --nproc-in 8 \
    --nproc-out 8 \
    --mark-dups \
    --output-stats {output.stats} \
    --output {output.pairsam} \
    {input}"""

# The pairtools split command is used to split the final .pairsam into two files: .sam (or .bam) and .pairs (.pairsam has two extra columns containing the alignments from which the Micro-C pair was extracted, these two columns are not included in .pairs files)
rule makePairsNSam:
    input:
        "results/rmPCRDup/{aliquot_barcode}.dedup.pairsam"
    output:
        sam = temp("results/makePairsNSam/{aliquot_barcode}.unsorted.sam"),
        pairs = "results/makePairsNSam/{aliquot_barcode}.mapped.pairs"
    threads: 16
    resources:
        mem_mb=32768
    conda:
        "micro-C"
    shell:"""
    pairtools split \
    --nproc-in 8 \
    --nproc-out 8 \
    --output-pairs {output.pairs} \
    --output-sam {output.sam} \
    {input}"""

rule finalBam:
    input:
        "results/makePairsNSam/{aliquot_barcode}.unsorted.sam"
    output:
        bam = "results/finalBam/{aliquot_barcode}.mapped.PT.bam",
        bai = "results/finalBam/{aliquot_barcode}.mapped.PT.bam.bai"
    threads: 32
    resources:
        mem_mb=65536
    shell:"""
    samtools sort -@16 \
    -T temp \
    -o {output.bam} \
    {input}

    samtools index {output.bam}"""

# Post-alignment QC

rule LibQC:
    input:
        "results/rmPCRDup/{aliquot_barcode}.stats.txt"
    output:
        "results/LibQC/{aliquot_barcode}.QC.txt"
    params:
        src = srcget_qc
    shell:"""
    {params.src} -p {input} > {output}"""

rule pairtoolsStats:
    input:
        "results/makePairsNSam/{aliquot_barcode}.mapped.pairs"
    output:
        "results/QC/{aliquot_barcode}.pairtools.stats.yaml"
    conda:
        "micro-C"
    shell:"""
    pairtools stats --yaml -o {output} {input}"""

rule samtoolsStats:
    input:
        "results/finalBam/{aliquot_barcode}.mapped.PT.bam"
    output:
        "results/QC/{aliquot_barcode}.samtools.stats.txt"
    shell:"""
    samtools stats {input} > {output}"""

rule multiQC:
    input:
        expand("results/QC/{aliquot_barcode}.pairtools.stats.yaml", aliquot_barcode=Sname),
        expand("results/QC/{aliquot_barcode}.samtools.stats.txt", aliquot_barcode=Sname),
        expand("results/rmPCRDup/{aliquot_barcode}.stats.txt", aliquot_barcode=Sname)
    output:
        "results/multiQC/multiqc_report.html"
    conda:
        "telomereC.py3.1"
    shell:"""
    multiqc results/ -o results/multiQC -f"""

# rule libComXQC:
#     input:
#         "results/finalBam/{aliquot_barcode}.mapped.PT.bam"
#     output:
#         "results/libComX/{aliquot_barcode}.libComX.preseq"
#     conda:
#         "micro-C"
#     shell:"""
#     preseq lc_extrap \
#     -bam \
#     -pe \
#     -extrap 2.1e9 \
#     -step 1e8 \
#     -seg_len 1000000000 \
#     -output {output} \
#     {input}"""

# Contact Maxrix and analysis

rule hicContacMatrix:
    input:
        "results/makePairsNSam/{aliquot_barcode}.mapped.pairs"
    output:
        "results/contacMatrix/{aliquot_barcode}.contact_map.hic"
    threads: 16
    resources:
        mem_mb=65536
    params:
        genome = ref_genome,
        src = srcJucier
    conda:
        "micro-C"
    shell:"""
    java -Xmx48000m  \
    -Djava.awt.headless=true \
    -jar {params.src} pre \
    --threads 16 \
    {input} \
    {output} \
    {params.genome}"""

rule mcoolContacMatrix:
    input:
        "results/makePairsNSam/{aliquot_barcode}.mapped.pairs"
    output:
        parisGZ = "results/makePairsNSam/{aliquot_barcode}.mapped.pairs.gz",
        cool = "results/coolerMatrix/{aliquot_barcode}.matrix_1kb.cool",
        mcool = "results/coolerMatrix/{aliquot_barcode}.matrix_1kb.mcool"
    params:
        genome = ref_genome,
        bin=1000
    threads: 32
    resources:
        mem_mb=65536
    conda:
        "micro-C"
    shell:"""
    bgzip -c {input} > {output.parisGZ}
    pairix {output.parisGZ}
    cooler cload pairix -p 16 {params.genome}:{params.bin} {output.parisGZ} {output.cool}
    cooler zoomify --balance -p 16 {output.cool}"""

rule bamCoverge:
    input:
        "results/finalBam/{aliquot_barcode}.mapped.PT.bam"
    output:
        "results/bamCoverage/{aliquot_barcode}.mapped.PT.bigwig"
    params:
        bin = 100,
        normalization = "RPKM"
    threads: 10
    conda:
        "telomereC.py3.1"
    shell:"""bamCoverage \
        -b {input} \
        -o {output} \
        -of bigwig \
        --binSize {params.bin} \
        --numberOfProcessors 10 \
        --ignoreDuplicates \
        --scaleFactor 1 \
        --normalizeUsing {params.normalization}
        """

rule parisDump:
    input:
        "results/makePairsNSam/{aliquot_barcode}.mapped.pairs"
    output:
        "results/makePairsNSam/{aliquot_barcode}.mapped.pairs.segdump.txt"
    params:
        conf = "config/circos.conf"
    shell:"""
        cat {input} | awk '{{print $1"\t"$2"\t"$3"\t"$3+1; print $1"\t"$4"\t"$5"\t"$5+1}}' | sed -n '/^#/!p' > {output}
        """

rule keepMultipleAlignPair2:
    input:
        "results/makePairsNSam/{aliquot_barcode}.mapped.pairs.segdump.txt"
    output:
        "results/makePairsNSam/{aliquot_barcode}.mapped.pairs.segdump.itinMA.txt"
    params:
        script="scripts/fixMAforCircos.py"
    shell:"""
        python {params.script} {input}
        """

# Feel free to adjust circos.template.conf for details pf ploting
#rule circoConfig2:
#    input:
#        "results/makePairsNSam/{aliquot_barcode}.mapped.pairs.segdump.itinMA.txt"
#    output:
#        new_conf = "results/plotCircos/{aliquot_barcode}.circos.config"
#    params:
#        karyotype = "/tgen_labs/barthel/software/miniforge3/envs/micro-C/data/karyotype/karyotype.human.chm13v2.txt",
#        tempConfig = "config/circos.template.conf",
#        script="scripts/makeCircosConf.sh"
#    shell:"""
#        cp {params.tempConfig} {output.new_conf}
#        sh {params.script} {params.karyotype} {input} {output.new_conf}
#    """

#rule plotCircos2:
#    input:
#        "results/plotCircos/{aliquot_barcode}.circos.config"
#    output:
#        "results/plotCircos/{aliquot_barcode}.circos.png",
#    conda:
#        "micro-C"
#    shell:"""
#        circos -conf {input} -noparanoid -outputfile {output}
#        """
