#!/bin/bash

module load hicpro
activate

cd /Hi-C/mapping
HiC-Pro -i /Hi-C/02.fasta/hdf_2_1 -o hdf_2_1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/hdf_2_2 -o hdf_2_2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/hdf_3_1 -o hdf_3_1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/hdf_3_2 -o hdf_3_2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/N_1_1 -o N_1_1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/N_1_2 -o N_1_2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/N_2_1 -o N_2_1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/N_2_2 -o N_2_2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/60M_L1 -o 60M_L1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/60M_L2 -o 60M_L2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/60M_S6_1 -o 60M_S6_1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/60M_S6_2 -o 60M_S6_2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/63F_L1 -o 63F_L1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/63F_L2 -o 63F_L2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/63F_S6_1 -o 63F_S6_1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/63F_S6_2 -o 63F_S6_2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/78M_L1 -o 78M_L1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/78M_L2 -o 78M_L2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/78M_S6_1 -o 78M_S6_1_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/78M_S6_2 -o 78M_S6_2_result -c /Hi-C/config-hicpro.txt -p 4
HiC-Pro -i /Hi-C/02.fasta/86F_L1 -o 86F_L1_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/86F_L2 -o 86F_L2_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/86F_S6_1 -o 86F_S6_1_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/86F_S6_2 -o 86F_S6_2_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/71F_L1 -o 71F_L1_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/71F_L2 -o 71F_L2_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/71F_S6_2 -o 71F_S6_2_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/71F_S6_2 -o 71F_S6_2_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/74F_L1 -o 74F_L1_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/74F_L2 -o 74F_L2_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/74F_S6_1 -o 74F_S6_1_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/74F_S6_2 -o 74F_S6_2_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/88M_L1 -o 88M_L1_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/88M_L2 -o 88M_L2_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/88M_S6_1 -o 88M_S6_1_result -c /Hi-C/config-hicpro.txt
HiC-Pro -i /Hi-C/02.fasta/88M_S6_2 -o 88M_S6_2_result -c /Hi-C/config-hicpro.txt
