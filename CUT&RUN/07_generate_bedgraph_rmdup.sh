#!/bin/bash


####https://yezhengstat.github.io/CUTTag_tutorial/#614_Visualization_of_peak_number,_peak_width,_peak_reproducibility_and_FRiPs
module load fastqc/0.11.8
module load trimmomatic/0.36
module load multiqc/1.11
module load bowtie2/2.3.5.1
module load samtools/1.9
module load bedtools/2.30.0

input="/deliv_Gorbunova_CUTandRUN_010324_raw"
projPath="/cutRun"
ref="//index_file/index/hg38_index"
chromSize="/index_file/hg38.chrom.sizes"



samplelist=(78M_L_Kap 78M_L_S6 78M_S6_Kap 78M_S6_S6 86F_L_Kap 86F_L_S6 86F_S6_Kap 86F_S6_S6 Positive Neo_S6 Neo)


cores=8


###############
for histName in "${samplelist[@]}"; do

#Spike-in calibration
#echo "$histName Normalized bedgraph"
#echo Job started at `date`

#seqDepthDouble=`samtools view -@ $cores -F 0x04 $projPath/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.sam | wc -l`
#seqDepth=$((seqDepthDouble/2))
#scale_factor=`echo "10000 / $seqDepth" | bc -l`
# echo "Scaling factor for $histName is: $scale_factor!"
#bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/rmDup/${histName}_bowtie2.rmDup.bed -g $chromSize &>$projPath/alignment/bedgraph/rmDup/${histName}_bowtie2.fragments.rmDup.normalized.bedgraph

bedtools genomecov -bg  -i $projPath/alignment/bed/rmDup/${histName}_bowtie2.rmDup.bed -g $chromSize &>$projPath/alignment/bedgraph/rmDup/${histName}_bowtie2.fragments.rmDup.bedgraph

done 
