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



outputname=(63F_L1 63F_L2 63F_S6_1 63F_S6_2 78MS6_1 78M_2 78M_L1 78M_S6_2 86F_L1 86F_L2 86F_S6_1 86F_S6_2 N_1_1 N_1_2 N_2_1 N_2_2 hdf_1_1 hdf_1_2 hdf_2_1 hdf_2_2 hdf_3_1 hdf_3_2)

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

