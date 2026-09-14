#!/bin/bash


####https://yezhengstat.github.io/CUTTag_tutorial/#614_Visualization_of_peak_number,_peak_width,_peak_reproducibility_and_FRiPs

module load bedtools/2.30.0
module load seacr/1.3
module load r/4.2.1/b1


input="/deliv_Gorbunova_CUTandRUN_010324_raw"
projPath="/cutRun"
ref="/index_file/index/hg38_index"
chromSize="/index_file/hg38.chrom.sizes"
seacr_path="/seacr/1.3/SEACR_1.3.sh"

mkdir -p $projPath/peakCalling/SEACR


samplelist=(78M_L_Kap 78M_L_S6 78M_S6_Kap 78M_S6_S6 86F_L_Kap 86F_L_S6 86F_S6_Kap 86F_S6_S6)

for histName in "${samplelist[@]}"; do 
    echo "$histName" 
    bash $seacr_path $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.eachSample.peaks
    echo "Completed: $histName"    
done
echo "All SEACR peak calling completed!"

