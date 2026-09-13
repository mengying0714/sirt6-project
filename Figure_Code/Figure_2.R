###########################################################################
#####    Code for Figure.2                                            #####
#####    Mengying Zhang (zhangmyharper@gmail.com)                     #####
###########################################################################

###########################################################################
#####    Code for Figure.2b                                          #####
###########################################################################

library(edgeR)
library(ggplot2)
library(ggrepel)
library(Rmisc)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)


HDF.ATAC.Peak.Count <- read.table("HDF.ATAC.counts.txt",header=T, sep="\t")
HDF.ATAC.Peak.Count <- subset(HDF.ATAC.Peak.Count, !(Chr == "chr19" & Start > 4174109 & End < 4182563))
dim(HDF.ATAC.Peak.Count)
# [1] 221950     16
HDF.ATAC.Peak.Count <- HDF.ATAC.Peak.Count[,c(1, 7:16)]
rownames(HDF.ATAC.Peak.Count) <- HDF.ATAC.Peak.Count$Peak
colnames(HDF.ATAC.Peak.Count) <- c("Peak", "S61F_Luc","S71F_Luc","S73M_Luc","S74F_Luc","S85F_Luc","S61F_SIRT6","S71F_SIRT6","S73M_SIRT6","S74F_SIRT6","S85F_SIRT6")
ATAC.Count <- HDF.ATAC.Peak.Count
dim(ATAC.Count)
# [1] 221953     11
ATAC.Count[is.na(ATAC.Count)] <- 0
ATAC.Count <- na.omit(ATAC.Count)
dim(ATAC.Count)
# [1] 221953     11
samples <- read.table("ATAC.HDF.samples.txt",header=T,sep="\t")
samples <- samples[,1:4]

y <- DGEList(counts=ATAC.Count[,2:11], group=samples$Group)
# colnames(y) <- samples$Sample
dim(y)
y <- calcNormFactors(y)
y$samples
### Data exploration
# plotMDS(y)
pdf("plotMDS.pdf", width=6, height=6)
plotMDS(y)
plotMDS(y, labels = y$samples$group)
dev.off()

y <- normLibSizes(y)
pdf("plotMDS after normLibSizes.pdf", width=6, height=6)
plotMDS(y)
plotMDS(y, labels = y$samples$group)
dev.off()

###  Estimating the dispersion
y <- estimateCommonDisp(y, verbose=TRUE)
# #Disp = 0.01503 , BCV = 0.1226
# ### Estimate gene-specific dispersions:
y <- estimateTagwiseDisp(y)
# plotBCV(y)

### Differential expression
#################### QL F-test ####################
### Design
## Paired
design <- model.matrix(~ subject + Group, data=samples)
y <- estimateDisp(y,design, robust=TRUE)
fit <- glmQLFit(y,design, robust=TRUE)
qlf <- glmQLFTest(fit, coef = 6)
topTags(qlf)

### Top edgeR_qlf results based QL F-test
#top <- topTags(qlf, n=nrow(y$counts))$table
top <- topTags(qlf, n=nrow(y$counts))$table; Count <- ATAC.Count[rownames(top),]; edgeR_qlf_Count <- cbind(top,Count)

## Add DEG infor 
edgeR_qlf_Count[which((edgeR_qlf_Count$FDR < 0.05) & (edgeR_qlf_Count$logFC > 0)),17]="1.UP"; edgeR_qlf_Count[which((edgeR_qlf_Count$FDR < 0.05) & (edgeR_qlf_Count$logFC < -0)),17]="0.Down"; colnames(edgeR_qlf_Count)[17] <- "ATAC.Diff"; edgeR_qlf_Count[is.na(edgeR_qlf_Count$ATAC.Diff),17]<- "2.Not"
write.table(edgeR_qlf_Count,"edgeR.HDF.ATAC.Peak.Count.Paired.txt",col.names=T,row.names=T,quote=F,sep="\t")

 
edgeR_qlf_Count <- edgeR_qlf_Count[4:221953,]
dim(edgeR_qlf_Count)
# [1] 221950     17

nrow(subset(edgeR_qlf_Count, FDR < 0.05 & logFC < 0));
nrow(subset(edgeR_qlf_Count, FDR < 0.05 & logFC > 0))
nrow(subset(edgeR_qlf_Count, FDR >= 0.05 & logFC < 0));
nrow(subset(edgeR_qlf_Count, FDR >= 0.05 & logFC > 0))
# > nrow(subset(edgeR_qlf_Count, FDR < 0.05 & logFC < 0));
# [1] 18452
# > nrow(subset(edgeR_qlf_Count, FDR < 0.05 & logFC > 0))
# [1] 7995
# > nrow(subset(edgeR_qlf_Count, FDR >= 0.05 & logFC < 0));
# [1] 94369
# > nrow(subset(edgeR_qlf_Count, FDR >= 0.05 & logFC > 0))
# [1] 101134

#### Volcano Plot
########## Volcano Plot of the difference between two knee groups

png(filename="Volcano Plot of the difference in ATAC peaks between two groups Paired.png", width=3.25,height=3.25,units="in",res=300)
ggplot(edgeR_qlf_Count, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ATAC.Diff), size = 0.2) +
  scale_color_manual(values = c("blue","red","black")) +
  ggtitle("Peak")+
  xlim(-5, 5)+
  ylim(0, 4.5)+
  geom_hline(yintercept=1.30103, linetype="dashed", color = "gray") +
  geom_vline(xintercept= 0, linetype="dashed", color = "gray") +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank())+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


###########################################################################
#####    Code for Figure.2c                                          #####
###########################################################################



edgeR_qlf_Count <- read.table("edgeR.HDF.ATAC.Peak.Count.Paired.txt",header=T, sep="\t")
HDF.ATAC.Peak.HMM <- read.table("HDF.ATAC.consensus_peaks.Closest.HMM.txt",header=F, sep="\t")
colnames(HDF.ATAC.Peak.HMM) <- c("Peak","HMM","Dis")

edgeR_qlf_Count.HMM <- inner_join(edgeR_qlf_Count, HDF.ATAC.Peak.HMM, by = "Peak")

write.table(edgeR_qlf_Count.HMM,"edgeR.HDF.ATAC.Peak.Count.Paired.HMM.txt",col.names=T,row.names=T,quote=F,sep="\t")

edgeR_qlf_Count.HMM.Up <- subset(edgeR_qlf_Count.HMM, ATAC.Diff == "1.UP")
edgeR_qlf_Count.HMM.Down <- subset(edgeR_qlf_Count.HMM, ATAC.Diff == "1.Down")
edgeR_qlf_Count.HMM.Other <- subset(edgeR_qlf_Count.HMM, ATAC.Diff == "2.Not")

edgeR_qlf_Count.HMM.Count <- edgeR_qlf_Count.HMM %>% count(ATAC.Diff, HMM)

write.table(edgeR_qlf_Count.HMM.Count,"edgeR.HDF.ATAC.Peak.Count.Paired.HMM.Count.txt",col.names=T,row.names=T,quote=F,sep="\t")

###########################################################################
#####    Code for Figure.2c                                          #####
###########################################################################

### Add CGI
edgeR_qlf_Count.HMM.TE <- read.table("edgeR.HDF.ATAC.Peak.Count.Paired.HMM.TE.txt",header=T, sep="\t")
HDF.ATAC.consensus_peaks.CGI <- read.table("HDF.ATAC.consensus_peaks.Closest.CGI.bed",header=F, sep="\t")
HDF.ATAC.consensus_peaks.CGI$V2 <- abs(HDF.ATAC.consensus_peaks.CGI$V2)
colnames(HDF.ATAC.consensus_peaks.CGI) <- c("Peak","Dis2CGI")
edgeR_qlf_Count.HMM.TE.CGI <- left_join(edgeR_qlf_Count.HMM.TE, HDF.ATAC.consensus_peaks.CGI, by = "Peak")

edgeR_qlf_Count.HMM.TE.CGI.CGI <- subset(edgeR_qlf_Count.HMM.TE.CGI, Dis2CGI == 0)
edgeR_qlf_Count.HMM.TE.CGI.CGIshore <- subset(edgeR_qlf_Count.HMM.TE.CGI, Dis2CGI != 0 & Dis2CGI > -2000 & Dis2CGI < 2000)
edgeR_qlf_Count.HMM.TE.CGI.Others <- subset(edgeR_qlf_Count.HMM.TE.CGI, Dis2CGI <= -2000 | Dis2CGI >= 2000)
edgeR_qlf_Count.HMM.TE.CGI.CGI$CGItype <- "CGI"
edgeR_qlf_Count.HMM.TE.CGI.CGIshore$CGItype <- "CGIshore"
edgeR_qlf_Count.HMM.TE.CGI.Others$CGItype <- "Others"
edgeR_qlf_Count.HMM.TE.CGI <- rbind(edgeR_qlf_Count.HMM.TE.CGI.CGI, edgeR_qlf_Count.HMM.TE.CGI.CGIshore, edgeR_qlf_Count.HMM.TE.CGI.Others)

write.table(edgeR_qlf_Count.HMM.TE.CGI,"edgeR.HDF.ATAC.Peak.Count.Paired.HMM.TE.CGI.txt",col.names=T,row.names=T,quote=F,sep="\t")

edgeR_qlf_Count.HMM.TE.CGI.Count <- edgeR_qlf_Count.HMM.TE.CGI %>% dplyr::count(ATAC.Diff, CGItype)

write.table(edgeR_qlf_Count.HMM.TE.CGI.Count,"edgeR.HDF.ATAC.Peak.Count.Paired.HMM.TE.CGI.Count.txt",col.names=T,row.names=T,quote=F,sep="\t")


###########################################################################
#####    Code for Figure.2f                                          #####
###########################################################################

load("edgeR.HDF.ATAC.Peak.Count.Paired.HMM.TE.CGI.TSS.DEG.RData")
TSS.DAR.HMM.TE.CGI.TSS.DEGs <- subset(edgeR_qlf_Count.HMM.TE.CGI.TSS.DEGs, Dis2TSS < 10000 & ATAC.Diff != "2.Not" & (SIRT6OE.type == "Up in SIRT6" | SIRT6OE.type == "Down in SIRT6"))
dim(TSS.DAR.HMM.TE.CGI.TSS.DEGs)

TSS.DAR.DEGs <- subset(TSS.DAR.HMM.TE.CGI.TSS.DEGs, Dis2TSS < 10000)
nrow(subset(TSS.DAR.DEGs, logFC > 0 & SIRT6OE.logFC > 0));
nrow(subset(TSS.DAR.DEGs, logFC > 0 & SIRT6OE.logFC < 0));
nrow(subset(TSS.DAR.DEGs, logFC < 0 & SIRT6OE.logFC > 0));
nrow(subset(TSS.DAR.DEGs, logFC < 0 & SIRT6OE.logFC < 0));
# > nrow(subset(TSS.DAR.DEGs, logFC > 0 & SIRT6OE.logFC > 0));
# [1] 494
# > nrow(subset(TSS.DAR.DEGs, logFC > 0 & SIRT6OE.logFC < 0));
# [1] 94
# > nrow(subset(TSS.DAR.DEGs, logFC < 0 & SIRT6OE.logFC > 0));
# [1] 211
# > nrow(subset(TSS.DAR.DEGs, logFC < 0 & SIRT6OE.logFC < 0));
# [1] 287

chromatin.organization.genes <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = "GO:0006325",
    columns = c("SYMBOL"),
    keytype = "GOALL"
)

unique(chromatin.organization.genes$SYMBOL)
TSS.DAR.DEGs.Chromatin <- subset(TSS.DAR.HMM.TE.CGI.TSS.DEGs, Dis2TSS < 10000 & Gene %in% chromatin.organization.genes$SYMBOL)
nrow(subset(TSS.DAR.DEGs.Chromatin, logFC > 0 & SIRT6OE.logFC > 0));
nrow(subset(TSS.DAR.DEGs.Chromatin, logFC > 0 & SIRT6OE.logFC < 0));
nrow(subset(TSS.DAR.DEGs.Chromatin, logFC < 0 & SIRT6OE.logFC > 0));
nrow(subset(TSS.DAR.DEGs.Chromatin, logFC < 0 & SIRT6OE.logFC < 0));
# > nrow(subset(TSS.DAR.DEGs.Chromatin, logFC > 0 & SIRT6OE.logFC > 0));
# [1] 46
# > nrow(subset(TSS.DAR.DEGs.Chromatin, logFC > 0 & SIRT6OE.logFC < 0));
# [1] 4
# > nrow(subset(TSS.DAR.DEGs.Chromatin, logFC < 0 & SIRT6OE.logFC > 0));
# [1] 10
# > nrow(subset(TSS.DAR.DEGs.Chromatin, logFC < 0 & SIRT6OE.logFC < 0));
# [1] 10

pdf(file='Compare difference between chromatin and gene changes after SIRT6OE.Dis2TSS10k only Chromatin genes with gene names.pdf', width=15, height=15)
ggplot(subset(TSS.DAR.DEGs.Chromatin, Gene != "SIRT6"), aes(x = logFC, y = SIRT6OE.logFC)) +
  geom_point(size = 1) +
  ggtitle("SIRT6OE")+
  xlim(-1, 2)+
  ylim(-1, 2)+
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept= 0, linetype="dashed", color = "gray") +
  geom_text(aes(label = Gene)) +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank())+
#   theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = "lm", se = FALSE) 
dev.off()
