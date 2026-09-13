#!/bin/bash

# Load necessary modules (adjust as needed for your system)
module load rsem/1.3.3

# Set paths
STAR_OUTPUT_DIR="/RNA_seq/fastq/STAR_New/"
RSEM_REF_DIR="/RNA_seq/fastq/rsem_index"
RSEM_OUTPUT_DIR="/RNA_seq/fastq/rsem"

# Create main output directory if it doesn't exist
mkdir -p $RSEM_OUTPUT_DIR

# Process all sample directories in the STAR output directory
for SAMPLE_DIR in $STAR_OUTPUT_DIR/*/
do

    SAMPLE_NAME=$(basename $SAMPLE_DIR)
    
    BAM_FILE="$SAMPLE_DIR/${SAMPLE_NAME}.Aligned.toTranscriptome.out.bam"
    if [ ! -f "$BAM_FILE" ]; then
        echo "Warning: $BAM_FILE not found. Skipping $SAMPLE_NAME"
        continue
    fi
    

    SAMPLE_OUTPUT_DIR="$RSEM_OUTPUT_DIR/$SAMPLE_NAME"
    mkdir -p $SAMPLE_OUTPUT_DIR
    
    # Run RSEM calculate expression
    rsem-calculate-expression --bam \
                              --no-bam-output \
                              -p 8 \
                              --paired-end \
                              $BAM_FILE \
                              $RSEM_REF_DIR/rsem_index \
                              $SAMPLE_OUTPUT_DIR/$SAMPLE_NAME
    
    echo "Processed $SAMPLE_NAME"
done

echo "All samples processed"