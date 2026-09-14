#!/bin/bash


####https://yezhengstat.github.io/CUTTag_tutorial/#614_Visualization_of_peak_number,_peak_width,_peak_reproducibility_and_FRiPs
module purge
module load bowtie2/2.3.5.1
module load samtools/1.11
module load deeptools/3.5.1

input="/deliv_Gorbunova_CUTandRUN_010324_raw"
projPath="/cutRun"
ref="/index_file/index/hg38_index"
chromSize="/index_file/hg38.chrom.sizes"



samplelist=(78M_L_Kap 78M_L_S6 78M_S6_Kap 78M_S6_S6 86F_L_Kap 86F_L_S6 86F_S6_Kap 86F_S6_S6 Positive)


cores=8
gsize=2913022398

###############
for histName in "${samplelist[@]}"; do

#####https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-5jyl8py98g2w/v2?step=14
samtools sort -o $projPath/alignment/bam/${histName}_bowtie2.mapped.sorted.bam $projPath/alignment/bam/${histName}_bowtie2.mapped.bam                                                     
samtools index $projPath/alignment/bam/${histName}_bowtie2.mapped.sorted.bam
#optional                                                                                                              


bamCoverage --normalizeUsing RPGC --effectiveGenomeSize ${gsize} -b $projPath/alignment/bam/${histName}_bowtie2.mapped.sorted.bam -o $projPath/alignment/bigwig/No_rmDUP/${histName}_mapped.sorted.normalized.bw

done
echo "All file convert completed!"