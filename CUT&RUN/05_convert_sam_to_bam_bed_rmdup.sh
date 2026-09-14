#!/bin/bash

####https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=12.2
###https://yezhengstat.github.io/CUTTag_tutorial/#614_Visualization_of_peak_number,_peak_width,_peak_reproducibility_and_FRiPs
module load fastqc/0.11.8
module load trimmomatic/0.36
module load multiqc/1.11
module load bowtie2/2.3.5.1
module load samtools/1.9
module load bedtools/2.30.0

input="/deliv_Gorbunova_CUTandRUN_010324_raw"
projPath="/cutRun"
ref="/index_file/index/hg38_index"
chromSize="/index_file/hg38.chrom.sizes"

samplelist=(78M_L_Kap 78M_L_S6 78M_S6_Kap 78M_S6_S6 86F_L_Kap 86F_L_S6 86F_S6_Kap 86F_S6_S6 Positive Neo_S6 Neo)


cores=8


for histName in "${samplelist[@]}"; do
seqDepthDouble=`samtools view -@ $cores -F 0x04 $projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/removeDuplicate/bowtie2_summary/${histName}_bowtie2.sorted.rmDup.seqDepth


##############
echo "$histName Convert Sam to Bam and Bed"
echo Job started at `date`
## Filter and keep the mapped read pairs
samtools view -bS -F 0x04 $projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam >$projPath/alignment/bam/rmDup/${histName}_bowtie2.rmDup.mapped.bam

## Convert into bed file format
bedtools bamtobed -i $projPath/alignment/bam/rmDup/${histName}_bowtie2.rmDup.mapped.bam -bedpe >$projPath/alignment/bed/rmDup/${histName}_bowtie2.rmDup.bed

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/rmDup/${histName}_bowtie2.rmDup.bed >$projPath/alignment/bed/rmDup/${histName}_bowtie2.rmDup.clean.bed

## Only extract the fragment related columns
cut -f 1,2,6 $projPath/alignment/bed/rmDup/${histName}_bowtie2.rmDup.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/rmDup/${histName}_bowtie2.rmDup.fragments.bed

done


multiqc /public/mzh103/MinseonCutTag/cutRun/alignment/removeDuplicate/*bowtie2.sorted.rmDup.sam /public/mzh103/MinseonCutTag/cutRun/alignment/removeDuplicate/