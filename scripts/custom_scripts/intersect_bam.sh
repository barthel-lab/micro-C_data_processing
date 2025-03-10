#!/bin/bash

# Dev
#inputs=(results/finalBam_filt/IMR90.filt.unmapped.PT.bam results/finalBam_filt/IMR90HiC.filt.unmapped.PT.bam)
#inputs=(results/finalBam_filt/BJHiC.filt.unmapped.PT.bam)
#peakBed="data/telomere-C_data/IMR90-capture.realn.mdup.MQ30.run_peaks.merge.bed"
#peakBed="data/telomere-C_data/BJP22-capture.realn.mdup.MQ30.run_peaks.merge.bed"
#peakBed="data/telomere-C_data/BJP14-capture.realn.mdup.MQ30.run_peaks.merge.bed"

# IMR90
inputs=(results/finalBam_filt/IMR90.filt.unmapped.PT.bam results/finalBam_filt/IMR90HiC.filt.unmapped.PT.bam)
peakBed="data/telomere-C_data/IMR90-capture.realn.mdup.MQ30.run_peaks.merge.bed"
ITS="data/telomere-C_data/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.TTAGGGs.simple.noTelomere.degenerated.bed"

for file in $(echo ${inputs[@]})
do name=$(basename -s .bam ${file[0]})
outbam="$name.intersect.bam"
outtxt="$name.intersect.IMR90.txt"
samtools view -b -q 30 $file | bedtools intersect -a stdin -b $peakBed > $outbam
echo "Number of Intersected pairs: $(samtools view $outbam |wc -l)" > $outtxt
echo "Number of Intersected read ID: $(samtools view $outbam | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Number of Intersected pairs in ITS: $(samtools view -b $outbam | bedtools intersect -a stdin -b $ITS | samtools view | wc -l)" >> $outtxt
echo "Number of Intersected read ID in ITS: $(samtools view -b $outbam | bedtools intersect -a stdin -b $ITS |samtools view | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Intersected read ID: " >> $outtxt
echo "$(samtools view $outbam | cut -f1 | sort -u)" >> $outtxt
done

# BJP22
inputs=(results/finalBam_filt/BJHiC.filt.unmapped.PT.bam)
peakBed="data/telomere-C_data/BJP22-capture.realn.mdup.MQ30.run_peaks.merge.bed"

for file in $(echo ${inputs[@]})
do name=$(basename -s .bam ${file[0]})
outbam="$name.intersect.bam"
outtxt="$name.intersect.BJP22.txt"
samtools view -b -q 30 $file | bedtools intersect -a stdin -b $peakBed > $outbam
echo "Number of Intersected pairs: $(samtools view $outbam |wc -l)" > $outtxt
echo "Number of Intersected read ID: $(samtools view $outbam | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Number of Intersected pairs in ITS: $(samtools view -b $outbam | bedtools intersect -a stdin -b $ITS | samtools view | wc -l)" >> $outtxt
echo "Number of Intersected read ID in ITS: $(samtools view -b $outbam | bedtools intersect -a stdin -b $ITS |samtools view | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Intersected read ID: " >> $outtxt
echo "$(samtools view $outbam | cut -f1 | sort -u)" >> $outtxt
done

# BJP14
inputs=(results/finalBam_filt/BJHiC.filt.unmapped.PT.bam)
peakBed="data/telomere-C_data/BJP14-capture.realn.mdup.MQ30.run_peaks.merge.bed"

for file in $(echo ${inputs[@]})
do name=$(basename -s .bam ${file[0]})
outbam="$name.intersect.bam"
outtxt="$name.intersect.BJP14.txt"
samtools view -b -q 30 $file | bedtools intersect -a stdin -b $peakBed > $outbam
echo "Number of Intersected pairs: $(samtools view $outbam |wc -l)" > $outtxt
echo "Number of Intersected read ID: $(samtools view $outbam | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Number of Intersected pairs in ITS: $(samtools view -b $outbam | bedtools intersect -a stdin -b $ITS | samtools view | wc -l)" >> $outtxt
echo "Number of Intersected read ID in ITS: $(samtools view -b $outbam | bedtools intersect -a stdin -b $ITS |samtools view | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Intersected read ID: " >> $outtxt
echo "$(samtools view $outbam | cut -f1 | sort -u)" >> $outtxt
done


