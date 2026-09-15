#!/bin/bash

####################################Deeploop###############################
###########                          hg38                  ################
###########                        2025.2.25               ################
###########                                                ################
###########################################################################
conda init
conda activate deeploop_py35

####prepare path and reference files
lib=/software/HiCorr/bin/Arima
# the path to lib of the Hi-C-data-preprocess 
#allValidPairs=/Hi-C/mapping/hdf_3_1_result/hic_results/data/hdf_3_1_result/hdf_3_1_result.allValidPairs
# from HiCPro
genome=hg38
# hg19/mm10
#enzyme=Arima
# restriction enzyme of the Hi-C experiment
#fragbed=/Hi-C/01.ref/annotion/filtered_frag_Arima.bed
# enzyme fragment bed "frag_1"..., provided in HiCorr reference files, HindIII, DPNII for mm10 and hg19 are provided, other type of reference files could be generated upon request
#outputname=hdf_3_1
#mkdir -p /Hi-C/Loop-deeploop/$outputname
# your output name or sample name
HiCorrPath=/software/HiCorr
# where you put "HiCorr" file
DeepLoopPath=/software/HiCorr/DeepLoop
# go to DeepLoop/
#DeepLoopBed=$DeepLoopPath/DeepLoop_models/ref/${genome}_${enzyme}_anchor_bed/


outputname=(71F_L1 71F_L2 71F_S6_1 71F_S6_2 74F_L1 74F_L2 74F_S6_1 74F_S6_2 88M_L1 88M_L2 88M_S6_1 88M_S6_2 HDF40F_1 HDF40F_2)

for outputname in "${outputname[@]}"; do

    echo "Processing sample: $outputname"
    
  #mkdir -p /Hi-C/Loop-deeploop/$outputname
  
  cd /Hi-C/Loop-deeploop/$outputname/


###run HICorr on cis and trans loop

#$HiCorrPath/HiCorr Arima frag_loop.$outputname.cis frag_loop.$outputname.trans $outputname $genome
#HiCorrOutputPath=`pwd`"/HiCorr_output/"
#DeepLoopOutputPath=`pwd`"/DeepLoop/"

#bash $bin/HindIII/HindIII.sh $ref/HindIII/ $bin/HindIII/ $2 $3 $4 $5  # <cis_loop_file> <trans_loop_file> <name_of_your_data> <reference_genome>
#bash $lib/HiCorr_Arima.sh  $ref/HindIII/ $lib/ frag_loop.$outputname.cis frag_loop.$outputname.trans $outputname $genome

  bash $lib/HiCorr_Arima.sh  /Hi-C/Loop-deeploop/Arima_ref/Arima $lib /Hi-C/Loop-deeploop/$outputname/frag_loop.$outputname.cis /Hi-C/Loop-deeploop/$outputname/frag_loop.$outputname.trans $outputname $genome
  
   wait
      echo "Completed sample: $outputname"
    echo "------------------------"
done

echo "All samples processed successfully"

