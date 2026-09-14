#!/usr/bin/env bash
# Reproducible paired-end ATAC-seq processing pipeline.
#
# Required environment variables:
#   PROJECT_DIR      Project directory containing raw/, bam/, bed/, bigwig/, peak/, and qc/
#   BOWTIE2_INDEX    Bowtie2 index prefix for the reference genome
#   REPEAT_BED       Sorted RepeatMasker BED file; column 4 must contain repName
#   CHROM_SIZES      Chromosome sizes file ordered like REPEAT_BED
#
# Optional environment variables:
#   SAMPLES          Space-separated sample names (default: study samples below)
#   PICARD_JAR       Path to picard.jar (default: /software/picard/3.4.0/picard.jar)
#   THREADS          Number of threads (default: 8)
#   GENOME_SIZE      MACS3 genome size (default: hs)
#   MODULES          Set to 0 if software is already available (default: 1)

#SBATCH --job-name=atac_seq
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --output=atac_seq_%j.log

set -euo pipefail

: "${PROJECT_DIR:?Set PROJECT_DIR before running this script.}"
: "${BOWTIE2_INDEX:?Set BOWTIE2_INDEX before running this script.}"
: "${REPEAT_BED:?Set REPEAT_BED before running this script.}"
: "${CHROM_SIZES:?Set CHROM_SIZES before running this script.}"

THREADS="${THREADS:-8}"
PICARD_JAR="${PICARD_JAR:-/software/picard/3.4.0/picard.jar}"
GENOME_SIZE="${GENOME_SIZE:-hs}"
MODULES="${MODULES:-1}"
DEFAULT_SAMPLES="61F_Luc 61F_SIRT6 71F_Luc 71F_SIRT6 73M_Luc 73M_SIRT6 74F_Luc 74F_SIRT6 85F_Luc 85F_SIRT6"
read -r -a samples <<< "${SAMPLES:-$DEFAULT_SAMPLES}"

RAW_DIR="${PROJECT_DIR}/raw"
QC_DIR="${PROJECT_DIR}/qc"
BAM_DIR="${PROJECT_DIR}/bam"
BED_DIR="${PROJECT_DIR}/bed"
BIGWIG_DIR="${PROJECT_DIR}/bigwig"
PEAK_DIR="${PROJECT_DIR}/peak"
COUNT_DIR="${PROJECT_DIR}/counts"
mkdir -p "$QC_DIR" "$BAM_DIR" "$BED_DIR" "$BIGWIG_DIR" "$PEAK_DIR" "$COUNT_DIR"

# Step 1: Load the software environment on a module-based HPC system.
if [[ "$MODULES" == "1" ]]; then
  module load fastqc/0.11.8 multiqc/1.21 bowtie2/2.3.5.1 samtools/1.21 picard/3.4.0 bedtools/2.30.0 deeptools/3.5.1 macs/3.0.3 subread/2.0.3
fi

process_sample() {
  local sample="$1"
  local read1="${RAW_DIR}/${sample}_R1.fastq.gz"
  local read2="${RAW_DIR}/${sample}_R2.fastq.gz"
  local sorted_bam="${BAM_DIR}/${sample}.mapped.sorted.bam"
  local rg_bam="${BAM_DIR}/${sample}.mapped.sorted.rg.bam"
  local final_bam="${BAM_DIR}/${sample}.rmdup.bam"
  local shifted_bam="${BAM_DIR}/${sample}.rmdup.shifted.bam"
  local sample_bed="${BED_DIR}/${sample}.rmdup.bed"
  local overlap_file="${COUNT_DIR}/${sample}.RepeatMasker.overlaps.tsv"
  local count_file="${COUNT_DIR}/${sample}.RepeatMasker.counts.tsv"
  local tmp_dir="${BAM_DIR}/${sample}.picard_tmp"

  [[ -f "$read1" && -f "$read2" ]] || { echo "ERROR: FASTQ files are missing for ${sample}." >&2; return 1; }

  # Step 2: Assess raw-read quality.
  fastqc --threads "$THREADS" --outdir "$QC_DIR" "$read1" "$read2"

  # Step 3: Align paired-end reads and retain primary mapped alignments.
  bowtie2 --threads "$THREADS" -x "$BOWTIE2_INDEX" -1 "$read1" -2 "$read2" 2> "${BAM_DIR}/${sample}.bowtie2.log" |
    samtools view -@ "$THREADS" -b -F 2052 - |
    samtools sort -@ "$THREADS" -o "$sorted_bam" -

  # Step 4: Add read-group metadata and remove PCR duplicates.
  java -Xmx10g -jar "$PICARD_JAR" AddOrReplaceReadGroups I="$sorted_bam" O="$rg_bam" \
    RGID="$sample" RGLB="${sample}_lib" RGPL=ILLUMINA RGPU=unit1 RGSM="$sample" VALIDATION_STRINGENCY=LENIENT
  rm -rf "$tmp_dir"
  mkdir -p "$tmp_dir"
  java -Xmx10g -jar "$PICARD_JAR" MarkDuplicates I="$rg_bam" O="$final_bam" \
    METRICS_FILE="${BAM_DIR}/${sample}.duplication_metrics.txt" REMOVE_DUPLICATES=true \
    CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT READ_NAME_REGEX=null TMP_DIR="$tmp_dir"

  # Step 5: Record alignment statistics.
  samtools flagstat "$final_bam" > "${BAM_DIR}/${sample}.rmdup.flagstat.txt"

  # Step 6: Apply the standard Tn5 insertion-site shift and index the result.
  alignmentSieve --bam "$final_bam" --ATACshift --numberOfProcessors "$THREADS" --outFile "$shifted_bam"
  samtools index "$shifted_bam"

  # Step 7: Convert alignments to BED and quantify RepeatMasker overlaps.
  bedtools bamtobed -i "$final_bam" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > "$sample_bed"
  bedtools intersect -sorted -a "$sample_bed" -b "$REPEAT_BED" -g "$CHROM_SIZES" -wa -wb |
    awk 'BEGIN{OFS="\t"}{print $1":"$2"-"$3,$7}' > "$overlap_file"
  { printf 'repName\tReadNum\n'; awk 'BEGIN{OFS="\t"}{count[$2]++} END{for (name in count) print name,count[name]}' "$overlap_file" | sort -k1,1; } > "$count_file"

  # Step 8: Generate an RPKM-normalized bigWig track from shifted fragments.
  bamCoverage --bam "$shifted_bam" --outFileName "${BIGWIG_DIR}/${sample}.rmdup.shifted.RPKM.bw" \
    --normalizeUsing RPKM --binSize 10 --numberOfProcessors "$THREADS"

  # Step 9: Call open-chromatin peaks with MACS3.
  macs3 callpeak -t "$shifted_bam" -f BAMPE -g "$GENOME_SIZE" -n "$sample" \
    --outdir "$PEAK_DIR" -q 0.01 --keep-dup all
}

for sample in "${samples[@]}"; do
  echo "[$(date)] Processing ${sample}"
  process_sample "$sample"
done

# Step 10: Summarize all FastQC reports.
multiqc "$QC_DIR" --outdir "$QC_DIR"

# Step 11: Merge sample peaks into a consensus peak set.
shopt -s nullglob
peak_files=("${PEAK_DIR}"/*_peaks.narrowPeak)
if (( ${#peak_files[@]} > 0 )); then
  cat "${peak_files[@]}" | sort -k1,1 -k2,2n | bedtools merge -i - > "${PEAK_DIR}/ATAC.consensus_peaks.bed"
  awk 'BEGIN{OFS="\t"}{print $1"."$2"."$3,$1,$2,$3,"."}' "${PEAK_DIR}/ATAC.consensus_peaks.bed" > "${PEAK_DIR}/ATAC.consensus_peaks.saf"

  # Step 12: Count paired-end ATAC fragments in consensus peaks.
  shifted_bams=("${BAM_DIR}"/*.rmdup.shifted.bam)
  featureCounts -a "${PEAK_DIR}/ATAC.consensus_peaks.saf" -F SAF -p --countReadPairs \
    -T "$THREADS" -o "${COUNT_DIR}/ATAC.consensus_peak_counts.txt" "${shifted_bams[@]}"

  # Optional Step 13: Count MeDIP fragments in ATAC peaks when MEDIP_BAM_DIR is provided.
  if [[ -n "${MEDIP_BAM_DIR:-}" ]]; then
    medip_bams=("${MEDIP_BAM_DIR}"/*.rmdup.bam)
    if (( ${#medip_bams[@]} > 0 )); then
      featureCounts -a "${PEAK_DIR}/ATAC.consensus_peaks.saf" -F SAF -p --countReadPairs \
        -T "$THREADS" -o "${COUNT_DIR}/ATAC_peaks.MeDIP_counts.txt" "${medip_bams[@]}"
    fi
  fi
fi

echo "[$(date)] ATAC-seq processing completed."
