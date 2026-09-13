###########################################################################
#####    Code for Figure.4                                            #####
#####    Mengying Zhang (zhangmyharper@gmail.com)                     #####
###########################################################################

###########################################################################
#####    Code for Figure.4b                                           #####
###########################################################################

#### Volcano Plot
########## Volcano Plot of the difference between two knee groups
library(ggplot2)
library(ggrepel)
library(Rmisc)
pdf(file='Volcano Plot of the difference between two groups Paired New.pdf', width=12, height=3)
p1 <- ggplot(edgeR_ATAC.LINE, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = ATAC.Diff), size = 1.5) +
  scale_color_manual(values = c("blue","black")) +
  ggtitle("LINE")+
  xlim(-0.4, 0.4)+
  ylim(0, 3)+
#   geom_hline(yintercept=1.30103, linetype="dashed", color = "gray") +
  geom_vline(xintercept= 0, linetype="dashed", color = "gray") +
  theme_bw(base_size = 16)+
  theme(panel.grid=element_blank())+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


###########################################################################
#####    Code for Figure.4c                                           #####
###########################################################################

computeMatrix scale-regions -S \
/Projects/CUTandTag_121223_result/alignment/bw/78M_L_ac_1_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/78M_L_ac_2_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/78M_L_me_1_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/78M_L_me_2_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/78M_S6_ac_1_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/78M_S6_ac_2_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/78M_S6_me_1_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/78M_S6_me_2_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/86F_L_ac_1_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/86F_L_ac_2_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/86F_L_me_1_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/86F_L_me_2_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/86F_S6_ac_1_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/86F_S6_ac_2_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/86F_S6_me_2_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/Luc_Neg_mapped.sorted.normalized.bw \
/Projects/CUTandTag_121223_result/alignment/bw/S6_Neg_mapped.sorted.normalized.bw \
-R SIRT6OE.Only.bed Control.Only.bed Both.bed \
--samplesLabel 78M_L_ac_1 78M_L_ac_2 78M_L_me_1 78M_L_me_2 78M_S6_ac_1 78M_S6_ac_2 78M_S6_me_1 78M_S6_me_2 86F_L_ac_1 86F_L_ac_2 86F_L_me_1 86F_L_me_2 86F_S6_ac_1 86F_S6_ac_2 86F_S6_me_2 Luc_Neg S6_Neg \
--beforeRegionStartLength 3000 \
--regionBodyLength 1000 \
--afterRegionStartLength 3000 \
--outFileName SIRT6Peak_CutTag_K9ac.K9me3_Matrix.3K.gz

plotHeatmap --colorList 'black, yellow, red' -m SIRT6Peak_CutTag_K9ac.K9me3_Matrix.3K_subset.gz -out SIRT6Peak_CutTag_K9ac.K9me3_Matrix.3K_subset_Order_by_RegionMean_BYR.png