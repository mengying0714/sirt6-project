###########################################################################
#####    Code for Figure.6                                            #####
#####    Mengying Zhang (zhangmyharper@gmail.com)                     #####
###########################################################################

###########################################################################
#####    Code for Figure.6c                                          #####
###########################################################################

library(dplyr)

# load("ATACSeq.Res.RData")

load("/Mouse/ATAC/Differential/ATACSeq.Res.RData")
load("/Mouse/MeDIP/Differential/MeDIP.Res.RData")

edgeR_MeQIP.Male.hSIRT6.VS.Luc$repName.1 <- rownames(edgeR_MeQIP.Male.hSIRT6.VS.Luc)

edgeR_ATAC.Male.hSIRT6.VS.Luc <- read.table("Male.hSIRT6.VS.Luc.edgeR_ATAC.Paired.txt",header=T,sep="\t")

edgeR_ATAC.res <- edgeR_ATAC.Male.hSIRT6.VS.Luc[,c("Family","SubFam","repName","repName.1","logFC","PValue","ATACSeq.Diff")]
colnames(edgeR_ATAC.res) <- c("Family","SubFam","repName","repName.1","ATAC.logFC","ATAC.PValue","ATAC.Diff")

edgeR_meDIP.res <- edgeR_MeQIP.Male.hSIRT6.VS.Luc[,c("repName.1","logFC","PValue","MeDIP.Diff")]
colnames(edgeR_meDIP.res) <- c("repName.1","MeDIP.logFC","MeDIP.PValue","MeDIP.Diff")

edgeR_ATAC.ATACSeq.res <- inner_join(edgeR_ATAC.res, edgeR_meDIP.res, by = "repName.1")
# edgeR_ATAC.ATACSeq.res <- inner_join(edgeR_ATAC.ATACSeq.res, TE.infor[,c(1,5)], by = "repName.1")
# edgeR_ATAC.ATACSeq.res$log10age <- log10(edgeR_ATAC.ATACSeq.res$age)
# edgeR_ATAC.ATACSeq.res$age.Billion <- edgeR_ATAC.ATACSeq.res$age / 1E9
# edgeR_ATAC.ATACSeq.res <- inner_join(edgeR_ATAC.ATACSeq.res, TE.Mean.Length, by = "repName.1")

edgeR_ATAC.ATACSeq.res.LINE <- subset(edgeR_ATAC.ATACSeq.res, Family == "LINE")
edgeR_ATAC.ATACSeq.res.LINE1 <- subset(edgeR_ATAC.ATACSeq.res, SubFam == "L1")

nrow(subset(edgeR_ATAC.ATACSeq.res.LINE1, MeDIP.logFC > 0 & ATAC.logFC > 0))
nrow(subset(edgeR_ATAC.ATACSeq.res.LINE1, MeDIP.logFC > 0 & ATAC.logFC < 0))
nrow(subset(edgeR_ATAC.ATACSeq.res.LINE1, MeDIP.logFC < 0 & ATAC.logFC > 0))
nrow(subset(edgeR_ATAC.ATACSeq.res.LINE1, MeDIP.logFC < 0 & ATAC.logFC < 0))
# > nrow(subset(edgeR_ATAC.ATACSeq.res.LINE1, MeDIP.logFC > 0 & ATAC.logFC > 0))
# [1] 71
# > nrow(subset(edgeR_ATAC.ATACSeq.res.LINE1, MeDIP.logFC > 0 & ATAC.logFC < 0))
# [1] 72
# > nrow(subset(edgeR_ATAC.ATACSeq.res.LINE1, MeDIP.logFC < 0 & ATAC.logFC > 0))
# [1] 1
# > nrow(subset(edgeR_ATAC.ATACSeq.res.LINE1, MeDIP.logFC < 0 & ATAC.logFC < 0))
# [1] 2
edgeR_MeQIP.Male.hSIRT6.VS.Luc.LINE <- subset(edgeR_MeQIP.Male.hSIRT6.VS.Luc, Family =="LINE")
pdf(file='Volcano Plot of the difference between edgeR_MeQIP.Male.hSIRT6.VS.Luc.LINE.With.Names.pdf', width=4.5, height=3)
ggplot(edgeR_MeQIP.Male.hSIRT6.VS.Luc.LINE, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = SubFam), size = 1.5) +
#   scale_color_manual(values = c("blue","red","black")) +
  ggtitle("LINE")+
  xlim(-1.3, 1.3)+
  ylim(0, 8)+
#   geom_hline(yintercept=1.30103, linetype="dashed", color = "gray") +
  geom_vline(xintercept= 0, linetype="dashed", color = "gray") +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank())+
#   theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()  


edgeR_ATAC.Male.hSIRT6.VS.Luc.LINE <- subset(edgeR_ATAC.Male.hSIRT6.VS.Luc, Family == "LINE")
pdf(file='Volcano Plot of the difference between edgeR_ATAC.Male.hSIRT6.VS.Luc.LINE.With.Names.pdf', width=4.5, height=3)
ggplot(edgeR_ATAC.Male.hSIRT6.VS.Luc.LINE, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = SubFam), size = 1.5) +
#   scale_color_manual(values = c("blue","red","black")) +
  ggtitle("LINE")+
  xlim(-0.55, 0.55)+
  ylim(0, 2)+
#   geom_hline(yintercept=1.30103, linetype="dashed", color = "gray") +
  geom_vline(xintercept= 0, linetype="dashed", color = "gray") +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank())+
#   theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()  
