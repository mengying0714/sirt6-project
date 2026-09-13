###########################################################################
#####    Code for Figure.1                                            #####
#####    Mengying Zhang (zhangmyharper@gmail.com)                     #####
###########################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggpubr)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Input data are not included; place the required files in the working directory before running this script.

###########################################################################
#####    Code for Figure.1c                                          #####
###########################################################################

select <- dplyr::select

# 1. Load and prepare TPM data
tpm_df.secondBatch <- read.table("merge.TPM.secondBatch.txt", 
                                  header=TRUE, sep="\t", check.names=FALSE)

cat("Number of columns:", ncol(tpm_df.secondBatch), "\n")
cat("Column names:", colnames(tpm_df.secondBatch), "\n")

tpm_df.secondBatch <- tpm_df.secondBatch %>%
    mutate(ensembl = str_extract(gene_id, "^ENSG\\d+"))

tpm_df.secondBatch$gene_symbol <- mapIds(org.Hs.eg.db,
                                          keys    = tpm_df.secondBatch$ensembl,
                                          keytype = "ENSEMBL",
                                          column  = "SYMBOL")

old_samples.secondBatch   <- c("61F_Luc_RNA","70F_Luc_RNA","71F_Luc_RNA","73M_Luc_RNA",
                                "74F_Luc_RNA","74M_Luc_RNA","79M_Luc_RNA","85F_Luc_RNA",
                                "86F_Luc_RNA","88M_Luc_RNA")
sirt6_samples.secondBatch <- c("61F_SIRT6_RNA","70F_SIRT6_RNA","71F_SIRT6_RNA","73M_SIRT6_RNA",
                                "74F_SIRT6_RNA","74M_SIRT6_RNA","79M_SIRT6_RNA","85F_SIRT6_RNA",
                                "86F_SIRT6_RNA","88M_SIRT6_RNA")

tpm_df1.secondBatch <- tpm_df.secondBatch %>% 
    select(all_of(c(old_samples.secondBatch, sirt6_samples.secondBatch)), gene_symbol)

expr_by_symbol.secondBatch <- tpm_df1.secondBatch %>%
    group_by(gene_symbol) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    filter(!is.na(gene_symbol))

cat("Number of genes:", nrow(expr_by_symbol.secondBatch), "\n")

# 2. Construct the expression matrix
expr_mat.secondBatch <- expr_by_symbol.secondBatch %>%
    column_to_rownames("gene_symbol") %>%
    dplyr::select(all_of(c(old_samples.secondBatch, sirt6_samples.secondBatch))) %>%
    as.matrix()

# 3. Filter low-expression genes
keep_genes.secondBatch <- 
    rowSums(expr_mat.secondBatch[, old_samples.secondBatch]   > 1) == length(old_samples.secondBatch) |
    rowSums(expr_mat.secondBatch[, sirt6_samples.secondBatch] > 1) == length(sirt6_samples.secondBatch)

expr_mat.secondBatch <- expr_mat.secondBatch[keep_genes.secondBatch, ]
cat("Number of genes retained after filtering:", nrow(expr_mat.secondBatch), "\n")

# 4. Log2 transformation
log2_mat.secondBatch <- log2(expr_mat.secondBatch + 1)

# 5. Construct the group and design matrix
donor.secondBatch <- gsub("_Luc_RNA|_SIRT6_RNA", "", 
                           c(old_samples.secondBatch, sirt6_samples.secondBatch))
donor.secondBatch <- factor(donor.secondBatch)

group.secondBatch <- factor(c(rep("Old",   length(old_samples.secondBatch)),
                               rep("SIRT6", length(sirt6_samples.secondBatch))),
                             levels = c("Old", "SIRT6"))

design.secondBatch <- model.matrix(~ donor.secondBatch + group.secondBatch)

# 6. Differential-expression analysis with limma
fit.secondBatch <- lmFit(log2_mat.secondBatch, design.secondBatch)
fit.secondBatch <- eBayes(fit.secondBatch)

result.secondBatch <- topTable(fit.secondBatch,
                                coef          = "group.secondBatchSIRT6",
                                number        = Inf,
                                sort.by       = "P",
                                adjust.method = "BH") %>%
    rownames_to_column("gene_symbol")

# 7. Optionally save the differential-expression results
write.csv(result.secondBatch, "DEG_SIRT6_vs_Old_secondBatch.csv", row.names = FALSE)
cat("Differential-expression analysis completed for", nrow(result.secondBatch), "genes\n")

# 8. Summarize significant genes
sig.secondBatch <- result.secondBatch %>%
    filter(adj.P.Val < 0.05) %>%
    mutate(type = ifelse(logFC > 0, "Up in SIRT6", "Down in SIRT6"))

cat("Significant genes (FDR < 0.05):", nrow(sig.secondBatch), "\n")
cat("  Upregulated in SIRT6:", sum(sig.secondBatch$logFC > 0), "\n")
cat("  Downregulated in SIRT6:", sum(sig.secondBatch$logFC < 0), "\n")

write.table(sig.secondBatch, "DEG_SIRT6_vs_Old_significant.txt",sep = "\t", quote = FALSE, row.names = FALSE)

# 9. Generate the volcano plot
volcano_df.secondBatch <- result.secondBatch %>%
    mutate(
        Significance = case_when(
            adj.P.Val < 0.05 & logFC >  0 ~ "Up in SIRT6",
            adj.P.Val < 0.05 & logFC <  0 ~ "Down in SIRT6",
            TRUE                           ~ "NS"
        )
    )

n_up.secondBatch   <- sum(volcano_df.secondBatch$Significance == "Up in SIRT6")
n_down.secondBatch <- sum(volcano_df.secondBatch$Significance == "Down in SIRT6")

options(repr.plot.width = 6, repr.plot.height = 4)

p <- ggplot(volcano_df.secondBatch, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = 0,            linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    annotate("text",
             x = max(volcano_df.secondBatch$logFC) * 0.9,
             y = max(-log10(volcano_df.secondBatch$adj.P.Val)) * 0.97,
             label = paste0("Up: ", n_up.secondBatch),
             color = "black", size = 5, hjust = 1) +
    annotate("text",
             x = min(volcano_df.secondBatch$logFC) * 0.9,
             y = max(-log10(volcano_df.secondBatch$adj.P.Val)) * 0.97,
             label = paste0("Down: ", n_down.secondBatch),
             color = "black", size = 5, hjust = 0) +
    scale_color_manual(values = c("Up in SIRT6"   = "#E74C3C",
                                  "Down in SIRT6" = "#3DBBD4",
                                  "NS"            = "grey70")) +
    theme_classic(base_size = 12) +
    theme(
        axis.text   = element_text(size = 15),
        axis.title  = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title  = element_text(size = 15)
    ) +
    labs(title = "SIRT6 vs Old — limma differential expression (Batch 2)",
         x     = "log2 Fold Change (SIRT6 / Old)",
         y     = "-log10(FDR)",
         color = NULL)

p
ggsave("Figure_1c_volcano_SIRT6_vs_Old_secondBatch.pdf", plot = p, width = 6, height = 4)


###########################################################################
#####    Code for Figure.1d                                          #####
###########################################################################

# Load previously saved GO enrichment results
go_up   <- read.table("GO_BP_OE_vs_Old_Up.txt",   header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
go_down <- read.table("GO_BP_OE_vs_Old_Down.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Plot enriched biological processes for upregulated genes
go_up_top <- go_up %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(log10P = -log10(p.adjust))

go_up_top$Description <- factor(go_up_top$Description,
                                 levels = go_up_top$Description[order(go_up_top$log10P)])

options(repr.plot.width = 6, repr.plot.height = 4)

p <- ggplot(go_up_top, aes(x = log10P, y = Description)) +
    geom_bar(stat = "identity", fill = "#E74C3C", width = 0.6) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "red", linewidth = 0.5) +
    theme_classic() +
    theme(axis.text.y   = element_text(size = 15),
          axis.text.x   = element_text(size = 15),
          axis.title.x  = element_text(size = 15),
          axis.title.y  = element_blank(),
          plot.title    = element_text(size = 15, hjust = 0.5)) +
    labs(x     = "-log10(FDR Q-Value)",
         title = "Up-regulated in Sirt6 OE\nGO Biological Process Enrichment")

p
ggsave("Figure_1d_GO_BP_OE_vs_Old_Up_bar.pdf", plot = p, width = 6, height = 4)

# Plot enriched biological processes for downregulated genes
go_down_top <- go_down %>%
    arrange(p.adjust) %>%
    head(15) %>%
    mutate(log10P = -log10(p.adjust))

go_down_top$Description <- factor(go_down_top$Description,
                                   levels = go_down_top$Description[order(go_down_top$log10P)])

options(repr.plot.width = 8, repr.plot.height = 4)

p <- ggplot(go_down_top, aes(x = log10P, y = Description)) +
    geom_bar(stat = "identity", fill = "#3d58a7", width = 0.6) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "red", linewidth = 0.5) +
    theme_classic() +
    theme(axis.text.y   = element_text(size = 13),
          axis.text.x   = element_text(size = 16),
          axis.title.x  = element_text(size = 16),
          axis.title.y  = element_blank(),
          plot.title    = element_text(size = 16, hjust = 0.5)) +
    labs(x     = "-log10(FDR Q-Value)",
         title = "Down-regulated in Sirt6 OE\nGO Biological Process Enrichment")

p
ggsave("Figure_1d_GO_BP_OE_vs_Old_Down_bar.pdf", plot = p, width = 8, height = 4)
     
###########################################################################
#####    Code for Figure.1e                                          #####
###########################################################################

# Reuse the TPM data loaded and annotated for Figure 1c
tpm_df2 <- tpm_df.secondBatch

options(repr.plot.width = 15, repr.plot.height = 4)

check_genes_query <- c("SUV39H1", "SUV39H2", "DNMT1", "DNMT3B")

plot_df2 <- tpm_df2 %>%
    filter(gene_symbol %in% check_genes_query) %>%
    dplyr::select(gene_symbol,
                  "61F_Luc_RNA","70F_Luc_RNA","71F_Luc_RNA","73M_Luc_RNA",
                  "74F_Luc_RNA","74M_Luc_RNA","79M_Luc_RNA","85F_Luc_RNA",
                  "86F_Luc_RNA","88M_Luc_RNA",
                  "61F_SIRT6_RNA","70F_SIRT6_RNA","71F_SIRT6_RNA","73M_SIRT6_RNA",
                  "74F_SIRT6_RNA","74M_SIRT6_RNA","79M_SIRT6_RNA","85F_SIRT6_RNA",
                  "86F_SIRT6_RNA","88M_SIRT6_RNA") %>%
    pivot_longer(cols = -gene_symbol,
                 names_to = "sample",
                 values_to = "TPM") %>%
    mutate(
        Group = case_when(
            grepl("Luc",   sample) ~ "Old",
            grepl("SIRT6", sample) ~ "SIRT6"
        ),
        Group = factor(Group, levels = c("Old", "SIRT6"))
    )

# Perform Wilcoxon tests and adjust for multiple comparisons
stat_res <- compare_means(TPM ~ Group,
                           data    = plot_df2,
                           group.by = "gene_symbol",
                           method  = "wilcox.test") %>%
    mutate(p.adj.BH = p.adjust(p, method = "BH"))


gene_order <- c("SUV39H1","SUV39H2","DNMT1","DNMT3B")
stat_res <- stat_res %>%
    mutate(gene_symbol = factor(gene_symbol, levels = gene_order),
           p.label = paste0("p.adj = ", formatC(p.adj.BH, format = "e", digits = 2)))

# Calculate label positions from the maximum TPM value of each gene
y_pos <- plot_df2 %>%
    group_by(gene_symbol) %>%
    summarise(y.position = max(TPM, na.rm = TRUE) * 1.1)

stat_res <- stat_res %>%
    left_join(y_pos, by = "gene_symbol")

plot_df2$gene_symbol <- factor(plot_df2$gene_symbol, levels = gene_order)

# Generate gene-expression boxplots
p <- ggplot(plot_df2, aes(x = Group, y = TPM, color = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.3, fill = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
    scale_color_manual(values = c("Old" = "#E74C3C", "SIRT6" = "#F5A623")) +
    facet_wrap(~ gene_symbol, scales = "free_y", nrow = 1) +
    geom_text(data = stat_res,
              aes(x = 1.5, y = y.position, label = p.label),
              inherit.aes = FALSE, size = 4) +
    theme_classic() +
    theme(legend.position  = "top",
          axis.title.x     = element_blank(),
          axis.text.x      = element_text(size = 15, angle = 45, hjust = 1),
          axis.text.y      = element_text(size = 15),
          axis.title.y     = element_text(size = 15),
          strip.text       = element_text(size = 15, face = "italic")) +
    labs(y = "TPM", color = "Group",
         title = "Batch 2: Old vs SIRT6")

p
ggsave("Figure_1e_boxplot_batch2_Old_vs_SIRT6.pdf", plot = p, width = 15, height = 4)
