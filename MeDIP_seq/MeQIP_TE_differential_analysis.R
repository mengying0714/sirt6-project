###############################################################################
# MeDIP-seq RepeatMasker RPKM and paired differential analysis
###############################################################################

library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(Rmisc)

data_dir <- Sys.getenv("MEDIP_DATA_DIR", unset=".")
output_dir <- Sys.getenv("MEDIP_OUTPUT_DIR", unset=".")
annotation_dir <- Sys.getenv("MEDIP_ANNOTATION_DIR", unset=".")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

###############################################################################
# Part 1: Read RepeatMasker read counts
###############################################################################

S61F_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.61F_Luc.Reads.Num.txt"), header=TRUE, sep="\t")
S71F_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.71F_Luc.Reads.Num.txt"), header=FALSE, sep=" ")
S73M_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.73M_Luc.Reads.Num.txt"), header=FALSE, sep=" ")
S74F_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.74F_Luc.Reads.Num.txt"), header=FALSE, sep=" ")
S85F_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.85F_Luc.Reads.Num.txt"), header=FALSE, sep=" ")
S61F_S6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.61F_S6.Reads.Num.txt"), header=FALSE, sep=" ")
S71F_S6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.71F_S6.Reads.Num.txt"), header=FALSE, sep=" ")
S73M_S6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.73M_S6.Reads.Num.txt"), header=FALSE, sep=" ")
S74F_S6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.74F_S6.Reads.Num.txt"), header=FALSE, sep=" ")
S85F_S6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.85F_S6.Reads.Num.txt"), header=FALSE, sep=" ")

# Assign consistent column names.
colnames(S61F_Luc) <- c("repName", "ReadNum")
colnames(S71F_Luc) <- c("repName", "ReadNum")
colnames(S73M_Luc) <- c("repName", "ReadNum")
colnames(S74F_Luc) <- c("repName", "ReadNum")
colnames(S85F_Luc) <- c("repName", "ReadNum")
colnames(S61F_S6) <- c("repName", "ReadNum")
colnames(S71F_S6) <- c("repName", "ReadNum")
colnames(S73M_S6) <- c("repName", "ReadNum")
colnames(S74F_S6) <- c("repName", "ReadNum")
colnames(S85F_S6) <- c("repName", "ReadNum")

# Replace missing values with zero.
S61F_Luc[is.na(S61F_Luc)] <- 0
S71F_Luc[is.na(S71F_Luc)] <- 0
S73M_Luc[is.na(S73M_Luc)] <- 0
S74F_Luc[is.na(S74F_Luc)] <- 0
S85F_Luc[is.na(S85F_Luc)] <- 0
S61F_S6[is.na(S61F_S6)] <- 0
S71F_S6[is.na(S71F_S6)] <- 0
S73M_S6[is.na(S73M_S6)] <- 0
S74F_S6[is.na(S74F_S6)] <- 0
S85F_S6[is.na(S85F_S6)] <- 0

###############################################################################
# Part 2: Add RepeatMasker length information
###############################################################################

RepName.Length <- read.table(file.path(annotation_dir, "RepeatMasker.hg38.Family.repName.merged.sorted.Total.Length.txt"), header=TRUE, sep="\t")

S61F_Luc.Length <- right_join(S61F_Luc, RepName.Length, by="repName")
S71F_Luc.Length <- right_join(S71F_Luc, RepName.Length, by="repName")
S73M_Luc.Length <- right_join(S73M_Luc, RepName.Length, by="repName")
S74F_Luc.Length <- right_join(S74F_Luc, RepName.Length, by="repName")
S85F_Luc.Length <- right_join(S85F_Luc, RepName.Length, by="repName")
S61F_S6.Length <- right_join(S61F_S6, RepName.Length, by="repName")
S71F_S6.Length <- right_join(S71F_S6, RepName.Length, by="repName")
S73M_S6.Length <- right_join(S73M_S6, RepName.Length, by="repName")
S74F_S6.Length <- right_join(S74F_S6, RepName.Length, by="repName")
S85F_S6.Length <- right_join(S85F_S6, RepName.Length, by="repName")

# Add the total number of duplicate-removed reads.
S61F_Luc.Length$TotalRead <- 90565248
S71F_Luc.Length$TotalRead <- 97546163
S73M_Luc.Length$TotalRead <- 85103859
S74F_Luc.Length$TotalRead <- 84607424
S85F_Luc.Length$TotalRead <- 88608546
S61F_S6.Length$TotalRead <- 76491031
S71F_S6.Length$TotalRead <- 92301779
S73M_S6.Length$TotalRead <- 89891755
S74F_S6.Length$TotalRead <- 83393166
S85F_S6.Length$TotalRead <- 88373325

# Convert read counts to numeric values.
S61F_Luc.Length$ReadNum <- as.numeric(S61F_Luc.Length$ReadNum)
S71F_Luc.Length$ReadNum <- as.numeric(S71F_Luc.Length$ReadNum)
S73M_Luc.Length$ReadNum <- as.numeric(S73M_Luc.Length$ReadNum)
S74F_Luc.Length$ReadNum <- as.numeric(S74F_Luc.Length$ReadNum)
S85F_Luc.Length$ReadNum <- as.numeric(S85F_Luc.Length$ReadNum)
S61F_S6.Length$ReadNum <- as.numeric(S61F_S6.Length$ReadNum)
S71F_S6.Length$ReadNum <- as.numeric(S71F_S6.Length$ReadNum)
S73M_S6.Length$ReadNum <- as.numeric(S73M_S6.Length$ReadNum)
S74F_S6.Length$ReadNum <- as.numeric(S74F_S6.Length$ReadNum)
S85F_S6.Length$ReadNum <- as.numeric(S85F_S6.Length$ReadNum)

###############################################################################
# Part 3: Calculate RPKM
###############################################################################

S61F_Luc.Length$S61F_Luc.RPKM <- S61F_Luc.Length$ReadNum / S61F_Luc.Length$TotalRead / S61F_Luc.Length$TotalLength * 1e9
S71F_Luc.Length$S71F_Luc.RPKM <- S71F_Luc.Length$ReadNum / S71F_Luc.Length$TotalRead / S71F_Luc.Length$TotalLength * 1e9
S73M_Luc.Length$S73M_Luc.RPKM <- S73M_Luc.Length$ReadNum / S73M_Luc.Length$TotalRead / S73M_Luc.Length$TotalLength * 1e9
S74F_Luc.Length$S74F_Luc.RPKM <- S74F_Luc.Length$ReadNum / S74F_Luc.Length$TotalRead / S74F_Luc.Length$TotalLength * 1e9
S85F_Luc.Length$S85F_Luc.RPKM <- S85F_Luc.Length$ReadNum / S85F_Luc.Length$TotalRead / S85F_Luc.Length$TotalLength * 1e9
S61F_S6.Length$S61F_S6.RPKM <- S61F_S6.Length$ReadNum / S61F_S6.Length$TotalRead / S61F_S6.Length$TotalLength * 1e9
S71F_S6.Length$S71F_S6.RPKM <- S71F_S6.Length$ReadNum / S71F_S6.Length$TotalRead / S71F_S6.Length$TotalLength * 1e9
S73M_S6.Length$S73M_S6.RPKM <- S73M_S6.Length$ReadNum / S73M_S6.Length$TotalRead / S73M_S6.Length$TotalLength * 1e9
S74F_S6.Length$S74F_S6.RPKM <- S74F_S6.Length$ReadNum / S74F_S6.Length$TotalRead / S74F_S6.Length$TotalLength * 1e9
S85F_S6.Length$S85F_S6.RPKM <- S85F_S6.Length$ReadNum / S85F_S6.Length$TotalRead / S85F_S6.Length$TotalLength * 1e9

###############################################################################
# Part 4: Combine RPKM and read-count tables
###############################################################################

MeDIP.RPKM <- left_join(S61F_Luc.Length[, c("repName", "S61F_Luc.RPKM")], S71F_Luc.Length[, c("repName", "S71F_Luc.RPKM")], by="repName")
MeDIP.RPKM <- left_join(MeDIP.RPKM, S73M_Luc.Length[, c("repName", "S73M_Luc.RPKM")], by="repName")
MeDIP.RPKM <- left_join(MeDIP.RPKM, S74F_Luc.Length[, c("repName", "S74F_Luc.RPKM")], by="repName")
MeDIP.RPKM <- left_join(MeDIP.RPKM, S85F_Luc.Length[, c("repName", "S85F_Luc.RPKM")], by="repName")
MeDIP.RPKM <- left_join(MeDIP.RPKM, S61F_S6.Length[, c("repName", "S61F_S6.RPKM")], by="repName")
MeDIP.RPKM <- left_join(MeDIP.RPKM, S71F_S6.Length[, c("repName", "S71F_S6.RPKM")], by="repName")
MeDIP.RPKM <- left_join(MeDIP.RPKM, S73M_S6.Length[, c("repName", "S73M_S6.RPKM")], by="repName")
MeDIP.RPKM <- left_join(MeDIP.RPKM, S74F_S6.Length[, c("repName", "S74F_S6.RPKM")], by="repName")
MeDIP.RPKM <- left_join(MeDIP.RPKM, S85F_S6.Length[, c("repName", "S85F_S6.RPKM")], by="repName")

MeDIP.Count <- left_join(S61F_Luc.Length[, c("repName", "ReadNum")], S71F_Luc.Length[, c("repName", "ReadNum")], by="repName")
MeDIP.Count <- left_join(MeDIP.Count, S73M_Luc.Length[, c("repName", "ReadNum")], by="repName")
MeDIP.Count <- left_join(MeDIP.Count, S74F_Luc.Length[, c("repName", "ReadNum")], by="repName")
MeDIP.Count <- left_join(MeDIP.Count, S85F_Luc.Length[, c("repName", "ReadNum")], by="repName")
MeDIP.Count <- left_join(MeDIP.Count, S61F_S6.Length[, c("repName", "ReadNum")], by="repName")
MeDIP.Count <- left_join(MeDIP.Count, S71F_S6.Length[, c("repName", "ReadNum")], by="repName")
MeDIP.Count <- left_join(MeDIP.Count, S73M_S6.Length[, c("repName", "ReadNum")], by="repName")
MeDIP.Count <- left_join(MeDIP.Count, S74F_S6.Length[, c("repName", "ReadNum")], by="repName")
MeDIP.Count <- left_join(MeDIP.Count, S85F_S6.Length[, c("repName", "ReadNum")], by="repName")

colnames(MeDIP.Count) <- c("repName", "S61F_Luc", "S71F_Luc", "S73M_Luc", "S74F_Luc", "S85F_Luc", "S61F_S6", "S71F_S6", "S73M_S6", "S74F_S6", "S85F_S6")
colnames(MeDIP.RPKM) <- c("repName", "S61F_Luc", "S71F_Luc", "S73M_Luc", "S74F_Luc", "S85F_Luc", "S61F_S6", "S71F_S6", "S73M_S6", "S74F_S6", "S85F_S6")

save(MeDIP.Count, MeDIP.RPKM, file=file.path(output_dir, "MeDIP.Res.RData"))

###############################################################################
# Part 5: Paired differential analysis using edgeR
###############################################################################

rownames(MeDIP.RPKM) <- MeDIP.RPKM$repName
rownames(MeDIP.Count) <- MeDIP.Count$repName
MeDIP.Count[is.na(MeDIP.Count)] <- 0
MeDIP.Count <- na.omit(MeDIP.Count)

samples <- read.table(file.path(data_dir, "MeDIP.HDF.samples.txt"), header=TRUE, sep="\t")
TE.infor <- read.table(file.path(annotation_dir, "TE.infor.txt"), header=TRUE, sep="\t")

y <- DGEList(counts=MeDIP.Count[, 2:11], group=samples$Group, lib.size=samples$lib.size)
y <- calcNormFactors(y)

# Generate MDS plots.
pdf(file.path(output_dir, "MeDIP_MDS.pdf"), width=6, height=6)
plotMDS(y)
plotMDS(y, labels=y$samples$group)
dev.off()

y <- normLibSizes(y)

pdf(file.path(output_dir, "MeDIP_MDS_after_normLibSizes.pdf"), width=6, height=6)
plotMDS(y)
plotMDS(y, labels=y$samples$group)
dev.off()

# Estimate dispersion.
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# Perform paired quasi-likelihood analysis.
design <- model.matrix(~ subject + Group, data=samples)
y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=6)
topTags(qlf)

# Extract all edgeR results.
top <- topTags(qlf, n=nrow(y$counts))$table
TE.infor1 <- TE.infor[rownames(top), ]
RPKM <- MeDIP.RPKM[rownames(top), ]
edgeR_qlf_RPKM <- cbind(TE.infor1, top, RPKM)

# Classify differential MeDIP signals using the original P-value threshold.
edgeR_qlf_RPKM$MeDIP.Diff <- "2.Not"
edgeR_qlf_RPKM$MeDIP.Diff[edgeR_qlf_RPKM$PValue < 0.05 & edgeR_qlf_RPKM$logFC > 0] <- "1.UP"
edgeR_qlf_RPKM$MeDIP.Diff[edgeR_qlf_RPKM$PValue < 0.05 & edgeR_qlf_RPKM$logFC < 0] <- "0.Down"

write.table(edgeR_qlf_RPKM, file.path(output_dir, "edgeR_MeDIP_Paired.txt"), col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

edgeR_MeQIP <- na.omit(edgeR_qlf_RPKM)
edgeR_MeQIP.LINE <- subset(edgeR_MeQIP, Family == "LINE")
edgeR_MeQIP.SINE <- subset(edgeR_MeQIP, Family == "SINE")
edgeR_MeQIP.LTR <- subset(edgeR_MeQIP, Family == "LTR")
edgeR_MeQIP.DNA <- subset(edgeR_MeQIP, Family == "DNA")

# Report the number of negative and positive log-fold changes.
nrow(subset(edgeR_MeQIP.LINE, logFC < 0)); nrow(subset(edgeR_MeQIP.LINE, logFC > 0))
nrow(subset(edgeR_MeQIP.SINE, logFC < 0)); nrow(subset(edgeR_MeQIP.SINE, logFC > 0))
nrow(subset(edgeR_MeQIP.LTR, logFC < 0)); nrow(subset(edgeR_MeQIP.LTR, logFC > 0))
nrow(subset(edgeR_MeQIP.DNA, logFC < 0)); nrow(subset(edgeR_MeQIP.DNA, logFC > 0))

save(edgeR_MeQIP, MeDIP.Count, MeDIP.RPKM, file=file.path(output_dir, "MeDIP.Res.RData"))

###############################################################################
# Part 6: Volcano plots by transposable-element class
###############################################################################

pdf(file.path(output_dir, "Volcano_Plot_Differential_MeDIP_Paired.pdf"), width=12, height=3)

p1 <- ggplot(edgeR_MeQIP.LINE, aes(x=logFC, y=-log10(PValue))) +
  geom_point(aes(color=MeDIP.Diff), size=1.5) +
  scale_color_manual(values=c("0.Down"="blue", "1.UP"="red", "2.Not"="black")) +
  ggtitle("LINE") + xlim(-0.2, 0.4) + ylim(0, 1.5) +
  geom_vline(xintercept=0, linetype="dashed", color="gray") +
  theme_bw(base_size=16) +
  theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(hjust=0.5))

p2 <- ggplot(edgeR_MeQIP.SINE, aes(x=logFC, y=-log10(PValue))) +
  geom_point(aes(color=MeDIP.Diff), size=1.5) +
  scale_color_manual(values=c("0.Down"="blue", "1.UP"="red", "2.Not"="black")) +
  ggtitle("SINE") + xlim(-0.2, 0.3) + ylim(0, 2) +
  geom_vline(xintercept=0, linetype="dashed", color="gray") +
  theme_bw(base_size=16) +
  theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(hjust=0.5))

p3 <- ggplot(edgeR_MeQIP.LTR, aes(x=logFC, y=-log10(PValue))) +
  geom_point(aes(color=MeDIP.Diff), size=1.5) +
  scale_color_manual(values=c("0.Down"="blue", "1.UP"="red", "2.Not"="black")) +
  ggtitle("LTR") + xlim(-0.4, 0.4) + ylim(0, 2) +
  geom_vline(xintercept=0, linetype="dashed", color="gray") +
  theme_bw(base_size=16) +
  theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(hjust=0.5))

p4 <- ggplot(edgeR_MeQIP.DNA, aes(x=logFC, y=-log10(PValue))) +
  geom_point(aes(color=MeDIP.Diff), size=1.5) +
  scale_color_manual(values=c("0.Down"="blue", "1.UP"="red", "2.Not"="black")) +
  ggtitle("DNA transposon") + xlim(-0.4, 0.4) + ylim(0, 2) +
  geom_vline(xintercept=0, linetype="dashed", color="gray") +
  theme_bw(base_size=16) +
  theme(panel.grid=element_blank(), legend.position="none", plot.title=element_text(hjust=0.5))

multiplot(p1, p2, p3, p4, cols=4)
dev.off()