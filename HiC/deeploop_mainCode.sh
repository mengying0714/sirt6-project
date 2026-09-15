#!/bin/bash

module load deeploop
activate

export CUDA_VISIBLE_DEVICES='-1'
export PYTHONPATH=$PYTHONPATH:/sfw2/deeploop/DeepLoop/

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_1/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_1/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_1/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_1/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "hdf_1 samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_2/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_2/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_2/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_2/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "hdf_2 samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_3/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_3/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_3/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/hdf_3/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "hdf_3 samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/HDF40F/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/HDF40F/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/HDF40F/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/HDF40F/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "HDF40F samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/63F/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/63F/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/63F/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/63F/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "63F samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/63F_S6/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/63F_S6/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/63F_S6/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/63F_S6/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "63F_S6 samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/78M/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/78M/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/78M/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/78M/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "78M samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/78M_S6/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/78M_S6/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/78M_S6/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/78M_S6/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "78M_S6 samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/86F/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/86F/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/86F/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/86F/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "86F samples processed successfully"


mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/86F_S6/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/86F_S6/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/86F_S6/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/86F_S6/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "86F_S6 samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/71F/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/71F/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/71F/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/71F/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "71F samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/71F_S6/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/71F_S6/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/71F_S6/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/71F_S6/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "71F_S6 samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/74F/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/74F/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/74F/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/74F/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "74F samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/74F_S6/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/74F_S6/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/74F_S6/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/74F_S6/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "74F_S6 samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/88M/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/88M/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/88M/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/88M/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "88M samples processed successfully"

mkdir -p /Hi-C/Loop-deeploop/mergeTechnicRep/88M_S6/deeploop_output
cd /Hi-C/Loop-deeploop/mergeTechnicRep/88M_S6/deeploop_output
for i in {1..22} X Y;do
  python /sfw2/deeploop/DeepLoop/prediction/predict_chromosome.py --full_matrix_dir /Hi-C/Loop-deeploop/mergeTechnicRep/88M_S6/HiCorr_output \
                                --input_name anchor_2_anchor.loop.chr${i} \
                                --h5_file /Hi-C/Loop-deeploop/Arima_ref/CPGZ_trained/100M.h5 \
                                --out_dir /Hi-C/Loop-deeploop/mergeTechnicRep/88M_S6/deeploop_output \
                                --anchor_dir /Hi-C/Loop-deeploop/Arima_ref/Arima/hg38_arima_anchor_bed/ \
                                --chromosome chr${i} --small_matrix_size 128  --step_size 128 --dummy 5 \
                                --val_cols obs exp
done
echo "88M_S6 samples processed successfully"












