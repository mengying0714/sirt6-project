###########################################################################
#####    Code for Figure.3                                            #####
#####    Mengying Zhang (zhangmyharper@gmail.com)                     #####
###########################################################################

###########################################################################
#####    Code for Figure.3a                                          #####
###########################################################################
library(VennDiagram)
library(RColorBrewer)
library(ggsci)

Control.Peaks <- read.table("HDF.SIRT6.CutRun.consensus_broadPeak.rmDup.q0.01.Control.bed",header=F, sep="\t")
SIRT6OE.Peaks <- read.table("HDF.SIRT6.CutRun.consensus_broadPeak.rmDup.q0.01.SIRT6OE.bed",header=F, sep="\t")
Control.Peaks$Peak <- paste(Control.Peaks$V1, Control.Peaks$V2, Control.Peaks$V3, sep = ".")
SIRT6OE.Peaks$Peak <- paste(SIRT6OE.Peaks$V1, SIRT6OE.Peaks$V2, SIRT6OE.Peaks$V3, sep = ".")
# myCol <- brewer.pal(2, "Pastel2") 
myCol <- pal_npg(palette = c("nrc"), alpha = 0.6)(2)  ## Nature color
# Chart
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")  ## suppress logging
venn.plot <- venn.diagram( # category.names = c("GE" , "ME", "GME"),
		x = list(Control = Control.Peaks$Peak, SIRT6OE = SIRT6OE.Peaks$Peak), 
        filename = NULL,                                    # Set NULL for pdf
        lwd = 2, lty = 'blank', fill = myCol,               # Circles     
        cex = .8, fontface = "bold", fontfamily = "sans",   # Numbers   
        cat.cex = 0.8, cat.fontface = "bold", cat.default.pos = "outer", cat.fontfamily = "sans", cat.pos = c(-27, 25), cat.dist = c(0.055, 0.06) # Set name
)
pdf(file="Venn Plot of Overlap of SIRT6 peaks.pdf", height = 1.5 , width = 1.7)
grid.draw(venn.plot)
dev.off()


###########################################################################
#####    Code for Figure.3b                                          #####
###########################################################################
library(dplyr)
library(edgeR)

cd /Human/CutRun/3.BroadPeak.to.ChromatinFeatures

awk 'NR>1{print $1"\t"$2"\t"$3}' SIRT6Peak_SIRT6_Binding_Matrix.3K.v2.sorted_regions.bed | sort -k1,1 -k2,2n > SIRT6Peak_SIRT6_Binding_Matrix.3K.v2.sorted_regions.new.bed

bedtools closest -a SIRT6Peak_SIRT6_Binding_Matrix.3K.v2.sorted_regions.new.bed -b /Users/su/Documents/Rochester/Lab/3.Projects/3.3.Epigenetics.Aging.HiC.Zhang/Human/ATAC/Peak_Diff/E055_15_coreMarks_hg38lift_dense.sorted.bed -D b -t first > SIRT6Peak_SIRT6_Binding_Matrix.3K.v2.sorted_regions.Closest.HMM.bed

cat SIRT6Peak_SIRT6_Binding_Matrix.3K.v2.sorted_regions.Closest.HMM.bed | awk '{print $1"."$2"."$3"\t"$7"\t"$13}' > SIRT6Peak_SIRT6_Binding_Matrix.3K.v2.sorted_regions.Closest.HMM.txt


HDF.SIRT6Peaks <- read.table("SIRT6Peak_SIRT6_Binding_Matrix.3K.v2.sorted_regions.txt",header=T, sep="\t")
colnames(HDF.SIRT6Peaks)[4] <- "Peak"
HDF.SIRT6Peaks$Peak <- paste(HDF.SIRT6Peaks$chrom, HDF.SIRT6Peaks$start, HDF.SIRT6Peaks$end,  sep = ".")
HDF.SIRT6Peaks <- HDF.SIRT6Peaks[,c(1:4, 11, 13)]

HDF.SIRT6.Peak.HMM <- read.table("SIRT6Peak_SIRT6_Binding_Matrix.3K.v2.sorted_regions.Closest.HMM.txt",header=F, sep="\t")
colnames(HDF.SIRT6.Peak.HMM) <- c("Peak","HMM","Dis")

HDF.SIRT6Peaks.HMM <- left_join(HDF.SIRT6Peaks, HDF.SIRT6.Peak.HMM, by = "Peak")
HDF.SIRT6Peaks.HMM[which(HDF.SIRT6Peaks.HMM$Dis > 0), "HMM"] <- "Other"

write.table(HDF.SIRT6Peaks.HMM,"edgeR.HDF.SIRT6.Peak.Count.Paired.HMM.txt",col.names=T,row.names=T,quote=F,sep="\t")

HDF.SIRT6Peaks.HMM.Count <- HDF.SIRT6Peaks.HMM %>% count(HMM, deepTools_group)

write.table(HDF.SIRT6Peaks.HMM.Count,"edgeR.HDF.SIRT6.Peak.Count.Paired.HMM.Count.txt",col.names=T,row.names=T,quote=F,sep="\t")


###########################################################################
#####    Code for Figure.3c                                           #####
###########################################################################
#####shell script
module load deeptools/3.5.1
cd /cutRun/Peak.Broad/computeMatrix
computeMatrix scale-regions -S \
/Projects/cutRun/alignment/bigwig/RPKM_rmDup/78M_L_S6_mapped.sorted.normalized.bw \
/Projects/cutRun/alignment/bigwig/RPKM_rmDup/86F_L_S6_mapped.sorted.normalized.bw \
/Projects/cutRun/alignment/bigwig/RPKM_rmDup/78M_S6_S6_mapped.sorted.normalized.bw \
/Projects/cutRun/alignment/bigwig/RPKM_rmDup/86F_S6_S6_mapped.sorted.normalized.bw \
/Projects/ChIP/ENCODE/NHDF.H3K4me1.ENCFF358WYW.bw \
/Projects/ChIP/ENCODE/NHDF.H3K4me2.ENCFF643ZVM.bw \
/Projects/ChIP/ENCODE/NHDF.H3K4me3.ENCFF953NOX.bw \
/Projects/ChIP/ENCODE/NHDF.H3K27ac.ENCFF754VWN.bw \
/Projects/ChIP/ENCODE/NHDF.H3K27me3.ENCFF785PLD.bw \
/Projects/ChIP/ENCODE/NHDF.H3K36me3.ENCFF837TBA.bw \
/Projects/ChIP/ENCODE/NHDF.H3K79me2.ENCFF129EQX.bw \
/Projects/ChIP/ENCODE/NHDF.H3K9ac.ENCFF859XFV.bw \
/Projects/ChIP/ENCODE/NHDF.H3K9me3.ENCFF993GEV.bw \
/Projects/ChIP/ENCODE/NHDF.H4K20me1.ENCFF290BBE.bw \
/Projects/ChIP/ENCODE/NHDF.CTCF.ENCFF007JHA.bw \
/Projects/ChIP/ENCODE/NHDF.H2A.Z.ENCFF841LDW.bw \
-R SIRT6OE.Only.bed Control.Only.bed Both.bed \
--samplesLabel 78M_L_S6 86F_L_S6 78M_S6_S6 86F_S6_S6 H3K4me1 H3K4me2 H3K4me3 H3K27ac H3K27me3 H3K36me3 H3K79me2 H3K9ac H3K9me3 H4K20me1 CTCF H2AZ \
--beforeRegionStartLength 3000 \
--regionBodyLength 1000 \
--afterRegionStartLength 3000 \
--outFileName SIRT6Peak_SIRT6_Binding_Matrix.ChIP.3K.v2.gz

#### Plot
plotHeatmap --colorList 'black, yellow, red' -m SIRT6Peak_SIRT6_Binding_Matrix.ChIP.3K.v2.gz --outFileSortedRegions SIRT6Peak_SIRT6_Binding_Matrix.ChIP.3K.v2.sorted_regions.bed -out SIRT6Peak_SIRT6_Binding_Matrix.ChIP.3K.v2_Order_by_RegionMean_BYR.png


###########################################################################
#####    Code for Figure.3d                                           #####
###########################################################################

#### Chi-squared test
read.table("edgeR.HDF.SIRT6.Peak.Count.Paired.Histone.Count.txt",header=T,sep="\t")-> ChIP_chi_square
ChIP_chi_square$Pvalue <- 1
ChIP_chi_square$estimate <- 0

for (i in 1:nrow(ChIP_chi_square)) {
    M <- as.table(rbind(c(ChIP_chi_square[i,2],ChIP_chi_square[i,4]), c(ChIP_chi_square[i,3],ChIP_chi_square[i,5])))
	M
	chisq.test(M)
	ChIP_chi_square[i, "Pvalue"] <- chisq.test(M)$p.value
	ChIP_chi_square[i, "estimate"] <- fisher.test(M)$estimate
}
write.table(ChIP_chi_square,"edgeR.HDF.SIRT6.Peak.Count.Paired.Histone.Count.ChIP_chi_square.txt",col.names=T,row.names=F,quote=F,sep="\t")
