#!/bin/bash


module load miniconda3/4.12.0
conda activate macs3


mkdir -p /cutRun/peakCalling/MACS3

macs3 callpeak -t /cutRun/alignment/bam/78M_S6_S6_bowtie2.mapped.bam \
    -f BAMPE -g hs \
    -n 78M_S6_S6 \
    --call-summits \
    --nomodel --extsize 150 \
    -q 0.05 \
    -B --SPMR \
    --keep-dup auto \
    --outdir /cutRun/peakCalling/MACS3 2> /cutRun/peakCalling/MACS3/78M_S6_S6_macs2.log
 
    
    
macs3 callpeak -t /cutRun/alignment/bam/78M_L_S6_bowtie2.mapped.bam \
    -f BAMPE -g hs \
    -n 78M_L_S6 \
    --call-summits \
    --nomodel --extsize 150 \
    -q 0.05 \
    -B --SPMR \
    --keep-dup auto \
    --outdir /cutRun/peakCalling/MACS3 2> /cutRun/peakCalling/MACS3/78M_L_S6_macs2.log
    
    
    
macs3 callpeak -t /cutRun/alignment/bam/86F_L_S6_bowtie2.mapped.bam \
    -f BAMPE -g hs \
    -n 86F_L_S6 \
    --call-summits \
    --nomodel --extsize 150 \
    -q 0.05 \
    -B --SPMR \
    --keep-dup auto \
    --outdir /cutRun/peakCalling/MACS3 2> /cutRun/peakCalling/MACS3/86F_L_S6_macs2.log
    
    
macs3 callpeak -t /cutRun/alignment/bam/86F_S6_S6_bowtie2.mapped.bam \
    -f BAMPE -g hs \
    -n 86F_S6_S6 \
    --call-summits \
    --nomodel --extsize 150 \
    -q 0.05 \
    -B --SPMR \
    --keep-dup auto \
    --outdir /cutRun/peakCalling/MACS3 2> /cutRun/peakCalling/MACS3/86F_S6_S6_macs2.log
    
macs3 callpeak -t /cutRun/alignment/bam/78M_L_Kap_bowtie2.mapped.bam \
    -f BAMPE -g hs \
    -n 78M_L_Kap \
    --call-summits \
    --nomodel --extsize 150 \
    -q 0.05 \
    -B --SPMR \
    --keep-dup auto \
    --outdir /cutRun/peakCalling/MACS3 2> /cutRun/peakCalling/MACS3/78M_L_Kap_macs2.log
    
    
macs3 callpeak -B -t /cutRun/alignment/bam/78M_S6_Kap_bowtie2.mapped.bam \
    -f BAMPE -g hs \
    -n 78M_S6_Kap \
    --call-summits \
    --nomodel --extsize 150 \
    -q 0.05 \
    -B --SPMR \
    --keep-dup auto \
    --outdir /cutRun/peakCalling/MACS3 2> /cutRun/peakCalling/MACS3/78M_S6_Kap_macs2.log
    
    
macs3 callpeak -B -t /cutRun/alignment/bam/86F_L_Kap_bowtie2.mapped.bam \
    -f BAMPE -g hs \
    -n 86F_L_Kap \
    --call-summits \
    --nomodel --extsize 150 \
    -q 0.05 \
    -B --SPMR \
    --keep-dup auto \
    --outdir /cutRun/peakCalling/MACS3 2> /cutRun/peakCalling/MACS3/86F_L_Kap_macs2.log
    
macs3 callpeak -B -t /cutRun/alignment/bam/86F_S6_Kap_bowtie2.mapped.bam \
    -f BAMPE -g hs \
    -n 86F_S6_Kap \
    --call-summits \
    --nomodel --extsize 150 \
    -q 0.05 \
    -B --SPMR \
    --keep-dup auto \
    --outdir /cutRun/peakCalling/MACS3 2> /cutRun/peakCalling/MACS3/86F_S6_Kap_macs2.log
    
    
    
