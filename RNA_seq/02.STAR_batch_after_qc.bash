#!/bin/bash


# Load necessary modules (adjust as needed for your system)
#module load STAR/2.7.3a
# Directory containing input fastq files (Trimmomatic output)
INPUT_DIR="/fastqc_trimmed_results012022"
# Directory for STAR output
OUTPUT_DIR="/RNA_seq/STAR_New/"
# STAR index directory
GENOME_DIR="/RNA_seq/STAR_index"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Process all R1 files
for R1_FILE in $INPUT_DIR/*_R1_val_1.fq.gz
do
   SAMPLE_NAME=$(basename $R1_FILE _R1_val_1.fq.gz)
   
   R2_FILE="${INPUT_DIR}/${SAMPLE_NAME}_R2_val_2.fq.gz"
   
   if [ ! -f "$R2_FILE" ]; then
       echo "Warning: $R2_FILE not found. Skipping $SAMPLE_NAME"
       continue
   fi
   
   echo "Processing sample: ${SAMPLE_NAME}"
   
   # Create sample-specific output directory
   SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$SAMPLE_NAME"
   mkdir -p $SAMPLE_OUTPUT_DIR
   
    # 运行STAR比对
    /software/STAR/STAR-2.7.9a/bin/Linux_x86_64/STAR \
        --runThreadN 16 \
        --genomeDir $GENOME_DIR \
        --readFilesIn $R1_FILE $R2_FILE \
        --readFilesCommand zcat \
        --outFileNamePrefix $SAMPLE_OUTPUT_DIR/$SAMPLE_NAME. \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN 16 \
        --quantMode TranscriptomeSAM GeneCounts
    
    echo "Finished processing $SAMPLE_NAME"
done

echo "All samples processed"

