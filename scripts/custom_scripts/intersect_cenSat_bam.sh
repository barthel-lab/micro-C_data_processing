#!/bin/bash

# Dev
#inputs=(results/finalBam_filt/IMR90.filt.unmapped.PT.bam results/finalBam_filt/IMR90HiC.filt.unmapped.PT.bam)
#inputs=(results/finalBam_filt/BJHiC.filt.unmapped.PT.bam)
#peakBed="data/telomere-C_data/IMR90-capture.realn.mdup.MQ30.run_peaks.merge.bed"
#peakBed="data/telomere-C_data/BJP22-capture.realn.mdup.MQ30.run_peaks.merge.bed"
#peakBed="data/telomere-C_data/BJP14-capture.realn.mdup.MQ30.run_peaks.merge.bed"
aSAT="/home/ychen/yc_data/chm13v2.0_censat_v2.1.3column.bed"

# IMR90
inputs=(results/finalBam_filt/IMR90.filt.unmapped.PT.bam results/finalBam_filt/IMR90HiC.filt.unmapped.PT.bam)

for file in $(echo ${inputs[@]})
do name=$(basename -s .bam ${file[0]})
outtxt="$name.intersect.IMR90.txt"

echo "Number of Intersected pairs in aSAT: $(samtools view -b $file | bedtools intersect -a stdin -b $aSAT | samtools view | wc -l)" > $outtxt
echo "Number of Intersected read ID in aSAT: $(samtools view -b $file | bedtools intersect -a stdin -b $aSAT |samtools view | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Intersected Sam: " >> $outtxt
echo "$(samtools view -b $file | bedtools intersect -a stdin -b $aSAT |samtools view)"  >> $outtxt

done

# BJ
inputs=(results/finalBam_filt/BJHiC.filt.unmapped.PT.bam)

for file in $(echo ${inputs[@]})
do name=$(basename -s .bam ${file[0]})
outtxt="$name.intersect.BJ.txt"

echo "Number of Intersected pairs in aSAT: $(samtools view -b $file | bedtools intersect -a stdin -b $aSAT | samtools view | wc -l)" > $outtxt
echo "Number of Intersected read ID in aSAT: $(samtools view -b $file | bedtools intersect -a stdin -b $aSAT |samtools view | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Intersected Sam: " >> $outtxt
echo "$(samtools view -b $file | bedtools intersect -a stdin -b $aSAT |samtools view)"  >> $outtxt

done


# For mapped
# IMR90
inputs=(results/finalBam_filt/IMR90.filt.mapped.PT.bam results/finalBam_filt/IMR90HiC.filt.mapped.PT.bam)

for file in $(echo ${inputs[@]})
do name=$(basename -s .bam ${file[0]})
outtxt="$name.intersect.IMR90.txt"

echo "Number of Intersected pairs in aSAT: $(samtools view -b $file | bedtools intersect -a stdin -b $aSAT | samtools view | wc -l)" > $outtxt
echo "Number of Intersected read ID in aSAT: $(samtools view -b $file | bedtools intersect -a stdin -b $aSAT |samtools view | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Intersected Sam: " >> $outtxt
echo "$(samtools view -b $file | bedtools intersect -a stdin -b $aSAT |samtools view)"  >> $outtxt

done

# BJ
inputs=(results/finalBam_filt/BJHiC.filt.mapped.PT.bam)

for file in $(echo ${inputs[@]})
do name=$(basename -s .bam ${file[0]})
outtxt="$name.intersect.BJ.txt"

echo "Number of Intersected pairs in aSAT: $(samtools view -b $file | bedtools intersect -a stdin -b $aSAT | samtools view | wc -l)" > $outtxt
echo "Number of Intersected read ID in aSAT: $(samtools view -b $file | bedtools intersect -a stdin -b $aSAT |samtools view | cut -f1 | sort -u | wc -l)" >> $outtxt

echo "Intersected Sam: " >> $outtxt
echo "$(samtools view -b $file | bedtools intersect -a stdin -b $aSAT |samtools view)"  >> $outtxt

done
