#!/bin/bash
#SBATCH --job-name=bowtie2
#SBATCH --output=bowtie2%j.out
#SBATCH --error=bowtie2%j.err
#SBATCH --time=2-0:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=preempt

####https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=12.2
module load fastqc/0.11.8
module load trimmomatic/0.36
module load multiqc/1.11
module load bowtie2/2.3.5.1
module load samtools/1.9


input="/deliv_Gorbunova_CUTandRUN_010324_raw"
projPath="/cutRun"
ref="/index_file/index/hg38_index"
#ADAPTER_FILE="/gpfs/fs1/sfw2/trimmomatic/0.36/adapters/TruSeq3-PE.fa"

#cat Neo1_S6_cut_R1.fastq.gz Neo2_S6_cut_R1.fastq.gz > Neo_S6_cut_R1.fastq.gz &
#cat Neo1_S6_cut_R2.fastq.gz Neo2_S6_cut_R2.fastq.gz > Neo_S6_cut_R2.fastq.gz &

cores=8
samplelist=(78M_L_Kap 78M_L_S6 78M_S6_Kap 78M_S6_S6 86F_L_Kap 86F_L_S6 86F_S6_Kap 86F_S6_S6 Positive)
## Build the bowtie2 reference genome index if needed:
## bowtie2-build path/to/hg38/fasta/hg38.fa /path/to/bowtie2Index/hg38
for hitsample in "${samplelist[@]}"; do

bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${ref} -1 ${projPath}/cutAdapter/${hitsample}_cut_R1.fastq.gz -2 ${projPath}/cutAdapter/${hitsample}_cut_R2.fastq.gz -S ${projPath}/alignment/sam/${hitsample}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${hitsample}_bowtie2.txt

done