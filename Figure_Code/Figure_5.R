###########################################################################
#####    Code for Figure.5                                            #####
#####    Mengying Zhang (zhangmyharper@gmail.com)                     #####
###########################################################################

library(tidyverse)
library(ggpubr)
library(limma)
library(preprocessCore)
library(pheatmap)

# Input data are not included; place the required files in the working directory before running this script.

###########################################################################
#####    Code for Figure.5c                                          #####
###########################################################################

df2 <- data.frame(sample = rep(1:6, 2),
                  group = c(rep("Old control", 6), rep("Old SIRT6", 6)),
                  value = c(15, 19, 23, 21, 21, 24, 14, 19, 21, 20, 19, 22))

df2$group <- factor(df2$group, levels = c("Old control", "Old SIRT6"))

p <- ggplot(df2, aes(x = group, y = value, group = sample)) +
  geom_line(color = "gray75", linewidth = 0.7) +
  geom_point(aes(color = group), size = 6) +
  scale_color_manual(values = c("Old control" = "#EF4B3E",
                                "Old SIRT6" = "#00A88F")) +
  stat_compare_means(aes(group = group), paired = TRUE, method = "t.test",
                     label = "p.format", label.x = 1.5, label.y = 26, size = 4) +
  scale_y_continuous(breaks = c(10, 15, 20, 25)) +
  coord_cartesian(ylim = c(10, 27)) +
  labs(x = NULL, y = "% of Inter-chromosome\ninteractions") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12))
p
ggsave("Figure_5c_interchromosomal_interactions.pdf", plot = p, width = 4, height = 4.5)
        
###########################################################################
#####    Code for Figure.5e                                          #####
###########################################################################

#########TE Family plot

old_gain <- read.table("CAT2_TE_family_enrichment_results.txt", header = TRUE, sep = "\t")
old_depleted <- read.table("CAT5_TE_family_enrichment_results.txt", header = TRUE, sep = "\t")

major_families <- c("Alu", "MIR", "SVA", "L1", "L2", "CR1", "RTE-BovB",
                    "ERVL-MaLR", "ERVL", "ERV1", "ERVK", "hAT-Charlie",
                    "TcMar-Tigger", "hAT-Tip100", "Simple_repeat", "Low_complexity")

# Merge old-gained and old-depleted enrichment results
combined_te <- old_gain %>%
  filter(family %in% major_families) %>%
  select(family, log2FC_old_gain = fold_change, p_adj_old_gain = p_adj) %>%
  left_join(old_depleted %>%
              filter(family %in% major_families) %>%
              select(family, log2FC_old_depleted = fold_change,
                     p_adj_old_depleted = p_adj),
            by = "family") %>%
  mutate(log2FC_old_gain = log2(log2FC_old_gain),
         log2FC_old_depleted = log2(log2FC_old_depleted),
         sig_old_gain = ifelse(p_adj_old_gain < 0.05, "sig", "ns"),
         sig_old_depleted = ifelse(p_adj_old_depleted < 0.05, "sig", "ns"),
         # Assign major TE classes
         te_class = case_when(
           family %in% c("Alu", "MIR", "SVA") ~ "SINE",
           family %in% c("L1", "L2", "CR1", "RTE-BovB") ~ "LINE",
           family %in% c("ERVL-MaLR", "ERVL", "ERV1", "ERVK") ~ "LTR",
           family %in% c("hAT-Charlie", "TcMar-Tigger", "hAT-Tip100") ~ "DNA",
           family %in% c("Simple_repeat", "Low_complexity") ~ "Other"),
         te_class = factor(te_class, levels = c("SINE", "LINE", "LTR", "DNA", "Other"))) %>%
  arrange(te_class, log2FC_old_gain)

# Check missing values
sum(is.na(combined_te))
print(combined_te %>% select(family, log2FC_old_gain, log2FC_old_depleted,
                             sig_old_gain, sig_old_depleted))

family_order <- combined_te$family

p <- ggplot(combined_te, aes(x = factor(family, levels = family_order))) +
  geom_segment(aes(y = log2FC_old_gain, yend = log2FC_old_depleted,
                   xend = factor(family, levels = family_order)),
               color = "grey70", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(aes(y = log2FC_old_gain, color = "Old gained",
                 alpha = sig_old_gain), size = 4, shape = 16) +
  geom_point(aes(y = log2FC_old_depleted, color = "Old depleted",
                 alpha = sig_old_depleted), size = 4, shape = 17) +
  facet_grid(. ~ te_class, scales = "free_x", space = "free_x") +
  scale_color_manual(values = c("Old gained" = "#D73027",
                                "Old depleted" = "#4575B4"), name = "") +
  scale_alpha_manual(values = c("sig" = 1, "ns" = 0.3),
                     name = "Adjusted P < 0.05",
                     labels = c("sig" = "Significant", "ns" = "NS")) +
  guides(alpha = guide_legend(override.aes = list(color = "black", size = 4)),
         color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  labs(y = expression(log[2]~"fold change vs random"), x = NULL,
       title = "TE family enrichment in old-gained and old-depleted loops") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = "right", panel.grid = element_blank(),
        plot.margin = margin(5, 40, 5, 5, "pt"))

p
ggsave("Figure_5e_TE_Old-Gained_vs_Old-Depleted_dumbbell.pdf", plot = p, width = 9, height = 4)

#########ChromHMM-state plot

# Load old-gained and old-depleted enrichment results
old_gain <- read.table("CAT2_chromHMMstatus_enrichment_results.txt", header = TRUE, sep = "\t")
old_depleted <- read.table("CAT5_chromHMMstatus_enrichment_results.txt", header = TRUE, sep = "\t")

# Merge enrichment results
combined <- old_gain %>%
  select(status, log2FC_old_gain = fold_change, p_adj_old_gain = p_adj) %>%
  left_join(old_depleted %>%
              select(status, log2FC_old_depleted = fold_change,
                     p_adj_old_depleted = p_adj),
            by = "status") %>%
  mutate(log2FC_old_gain = log2(log2FC_old_gain),
         log2FC_old_depleted = log2(log2FC_old_depleted),
         sig_old_gain = ifelse(p_adj_old_gain < 0.05, "sig", "ns"),
         sig_old_depleted = ifelse(p_adj_old_depleted < 0.05, "sig", "ns"),
         # Assign ChromHMM classes
         hmm_class = case_when(
           status %in% c("1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk") ~ "Promoter",
           status %in% c("3_TxFlnk", "4_Tx", "5_TxWk") ~ "Transcription",
           status %in% c("6_EnhG", "7_Enh", "12_EnhBiv") ~ "Enhancer",
           status %in% c("9_Het", "8_ZNF/Rpts") ~ "Heterochromatin",
           status %in% c("13_ReprPC", "14_ReprPCWk") ~ "Polycomb",
           status == "15_Quies" ~ "Quiescent"),
         hmm_class = factor(hmm_class, levels = c("Promoter", "Transcription",
                                                  "Enhancer", "Polycomb",
                                                  "Heterochromatin", "Quiescent"))) %>%
  arrange(hmm_class, log2FC_old_gain)

# Set the ChromHMM-state order
status_order <- combined$status
options(repr.plot.width = 9, repr.plot.height = 4)

p <- ggplot(combined, aes(x = factor(status, levels = status_order))) +
  geom_segment(aes(y = log2FC_old_gain, yend = log2FC_old_depleted,
                   xend = factor(status, levels = status_order)),
               color = "grey70", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(aes(y = log2FC_old_gain, color = "Old gained",
                 alpha = sig_old_gain), size = 4, shape = 16) +
  geom_point(aes(y = log2FC_old_depleted, color = "Old depleted",
                 alpha = sig_old_depleted), size = 4, shape = 17) +
  facet_grid(. ~ hmm_class, scales = "free_x", space = "free_x") +
  scale_color_manual(values = c("Old gained" = "#D73027",
                                "Old depleted" = "#4575B4"), name = "") +
  scale_alpha_manual(values = c("sig" = 1, "ns" = 0.3),
                     name = "Adjusted P < 0.05",
                     labels = c("sig" = "Significant", "ns" = "NS")) +
  guides(alpha = guide_legend(override.aes = list(color = "black", size = 4)),
         color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  labs(y = expression(log[2]~"fold change vs random"), x = NULL,
       title = "ChromHMM enrichment in old-gained and old-depleted loops") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(5, 40, 5, 5, "pt"))

p
ggsave("Figure_5e_ChromHMM_Old-Gained_vs_Old-Depleted_dumbbell.pdf",plot = p, width = 9, height = 4)


###########################################################################
#####    Code for Figure.5g                                          #####
###########################################################################

# Load Old and SIRT6-OE loop matrices
oldE <- read.table("merged_oldGroup.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
oeE <- read.table("merged_OEGroup.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

old_samples <- c("63F", "71F", "74F", "78M", "86F", "88M")
oe_samples <- c("63F_S6", "71F_S6", "74F_S6", "78M_S6", "86F_S6", "88M_S6")

# Generate loop IDs and retain sample columns
oldE$combined_id <- apply(oldE[, 1:6], 1, paste, collapse = "_")
oeE$combined_id <- apply(oeE[, 1:6], 1, paste, collapse = "_")

result_old <- oldE[, 7:ncol(oldE)]
result_oe <- oeE[, 7:ncol(oeE)]
colnames(result_old) <- c(old_samples, "combined_id")
colnames(result_oe) <- c(oe_samples, "combined_id")

# Merge Old and OE matrices
merged_data <- full_join(result_old, result_oe, by = "combined_id")
merged_data1 <- merged_data[, c(old_samples, oe_samples)]
rownames(merged_data1) <- merged_data$combined_id
merged_data1[is.na(merged_data1)] <- 0

mat <- data.matrix(merged_data1)
cat("Samples included in the analysis:\n")
print(colnames(mat))

# Quantile normalization followed by log transformation
mat_qnorm <- normalize.quantiles(mat)
rownames(mat_qnorm) <- rownames(mat)
colnames(mat_qnorm) <- colnames(mat)
mat_norm <- log1p(mat_qnorm)

# Construct the design matrix for Old versus SIRT6-OE
group <- factor(c(rep("Old", length(old_samples)), rep("OE", length(oe_samples))),levels = c("Old", "OE"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Differential loop analysis using limma
fit <- lmFit(mat_norm, design)
contrast_matrix <- makeContrasts(OE_vs_Old = OE - Old, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))
res_OE_vs_Old <- topTable(fit2, coef = "OE_vs_Old", number = Inf, adjust.method = "BH")

# Identify increased and decreased loops
sig_OE_up <- res_OE_vs_Old[res_OE_vs_Old$P.Value < 0.001 & res_OE_vs_Old$logFC > 0, ]
sig_OE_down <- res_OE_vs_Old[res_OE_vs_Old$P.Value < 0.001 & res_OE_vs_Old$logFC < 0, ]

cat("Number of increased loops in SIRT6-OE:", nrow(sig_OE_up), "\n")
cat("Number of decreased loops in SIRT6-OE:", nrow(sig_OE_down), "\n")

write.table(sig_OE_up, "sig_OE_vs_Old.OEup.limma.significant.v4.txt",sep = "\t", quote = FALSE, row.names = TRUE)
write.table(sig_OE_down, "sig_OE_vs_Old.OEdown.limma.significant.v4.txt",sep = "\t", quote = FALSE, row.names = TRUE)

# Prepare significant loops for heatmap visualization
sig_ids <- rownames(res_OE_vs_Old)[res_OE_vs_Old$P.Value < 0.001]
heatmap_data <- mat_norm[sig_ids, c(old_samples, oe_samples), drop = FALSE]
heatmap_data_z <- t(scale(t(heatmap_data)))
heatmap_data_z <- heatmap_data_z[complete.cases(heatmap_data_z), , drop = FALSE]

annotation_col <- data.frame(Group = factor(group, levels = c("Old", "OE")))
rownames(annotation_col) <- colnames(heatmap_data_z)
ann_colors <- list(Group = c("Old" = "#4575B4", "OE" = "#D73027"))

pheatmap(heatmap_data_z, cluster_rows = TRUE, cluster_cols = FALSE,annotation_col = annotation_col, annotation_colors = ann_colors,show_rownames = FALSE, color = colorRampPalette(c("blue", "white", "red"))(100),fontsize = 12, border_color = NA,filename = "Figure_5g_Old_vs_SIRT6OE_significant_loops_heatmap.pdf",width = 6, height = 8)


###########################################################################
#####    Code for Figure.5h                                          #####
###########################################################################

# Load GREAT biological-process enrichment results
great_res <- read.table("sig_OE_vs_Old.OEDown.selected.BP.txt",header = TRUE, sep = "\t")
# Calculate enrichment significance and order biological processes
great_res$log10FDR <- -log10(great_res$BinomFdrQ)
great_res$Desc <- factor(great_res$Desc,levels = great_res$Desc[order(great_res$log10FDR)])
options(repr.plot.width = 11, repr.plot.height = 6)
p <- ggplot(great_res, aes(x = log10FDR, y = Desc)) +
  geom_col(fill = "#2E86C1", width = 0.6) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             color = "red", linewidth = 0.5) +
  labs(x = expression(-log[10]("FDR Q-value")), y = NULL,
       title = "Biological process enrichment of loops decreased in SIRT6 OE") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5))
p
ggsave("Figure_5h_SIRT6OE_decreased_loops_GO_BP_enrichment.pdf", plot = p, width = 11, height = 6)


###########################################################################
#####    Code for Figure.5i                                          #####
###########################################################################

# Load GREAT biological-process enrichment results
great_res <- read.table("sig_OE_vs_Old.OEup.selected.BP.txt", header = TRUE, sep = "\t")

# Calculate enrichment significance and order biological processes
great_res$log10FDR <- -log10(great_res$BinomFdrQ)
great_res$Desc <- factor(great_res$Desc,levels = great_res$Desc[order(great_res$log10FDR)])

options(repr.plot.width = 9, repr.plot.height = 6)

p <- ggplot(great_res, aes(x = log10FDR, y = Desc)) +
  geom_col(fill = "#2E86C1", width = 0.6) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             color = "red", linewidth = 0.5) +
  labs(x = expression(-log[10]("FDR Q-value")), y = NULL,
       title = "Biological process enrichment of loops increased in SIRT6 OE") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5))
p
ggsave("Figure_5i_SIRT6OE_increased_loops_GO_BP_enrichment.pdf", plot = p, width = 9, height = 6)
