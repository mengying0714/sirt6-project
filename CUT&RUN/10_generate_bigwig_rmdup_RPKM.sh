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



samplelist=(78M_L_S6 78M_S6_S6 86F_L_S6 86F_S6_S6)




###############
for histName in "${samplelist[@]}"; do

                                                                                                          


bamCoverage --normalizeUsing RPKM -b $projPath/alignment/bam/rmDup/${histName}_bowtie2.rmDup.mapped.sorted.bam -o $projPath/alignment/bigwig/RPKM_rmDup/${histName}_mapped.sorted.normalized.bw

done
echo "All file convert completed!"