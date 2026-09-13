#!/bin/bash


module load cutadapt/b2
module load trimgalore/0.6.2

INPUT_DIR=/deliv_RNA_012022_raw
OUTPUT_DIR=/fastqc_trimmed_results012022
mkdir -p ${OUTPUT_DIR}

for R1_FILE in ${INPUT_DIR}/*_R1.fastq.gz; do
    BASE_NAME=$(basename ${R1_FILE} _R1.fastq.gz)
    R2_FILE="${R1_FILE%_R1.fastq.gz}_R2.fastq.gz"
    
    if [ ! -f "$R2_FILE" ]; then
        echo "Warning: $R2_FILE not found. Skipping $BASE_NAME"
        continue
    fi
    
    echo "Processing sample: ${BASE_NAME}"
    
    trim_galore \
        -q 25 \
        --phred33 \
        --length 35 \
        --stringency 3 \
        --nextera \
        --paired \
        --clip_R1 15 \
        --clip_R2 15 \
        -o ${OUTPUT_DIR} \
        ${R1_FILE} \
        ${R2_FILE}
    
    echo "Finished processing sample: ${BASE_NAME}"
    echo "Output saved to: ${OUTPUT_DIR}"
done