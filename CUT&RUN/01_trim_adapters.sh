#!/bin/bash


####https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=12.2
module load fastqc/0.11.8
module load trimmomatic/0.36
module load cutadapt/b2

input="/deliv_Gorbunova_CUTandRUN_010324_raw"
projPath="/cutRun"
ref="/index_file/index/hg38_index"
#ADAPTER_FILE="/trimmomatic/0.36/adapters/TruSeq3-PE.fa"

samplelist=(78M_L_Kap 78M_L_S6 78M_S6_Kap 78M_S6_S6 86F_L_Kap 86F_L_S6 86F_S6_Kap 86F_S6_S6 Positive)

for hitsample in "${samplelist[@]}"; do

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 20 -o ${projPath}/cutAdapter/${hitsample}_cut_R1.fastq.gz -p ${projPath}/cutAdapter/${hitsample}_cut_R2.fastq.gz ${input}/${hitsample}_R1.fastq.gz ${input}/${hitsample}_R2.fastq.gz
done

echo "cutAdapter finished!"