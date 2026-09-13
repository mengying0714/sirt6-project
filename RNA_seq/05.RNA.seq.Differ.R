###############################################################################
# Prepare and analyze RNA-seq TPM data
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggpubr)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Run this script from the directory containing the sample folders.
base_dir <- "."
output_dir <- "/TPM"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

donors <- c("61F", "70F", "71F", "73M", "74F", "74M", "79M", "85F", "86F", "88M")
samples_old <- paste0(donors, "_Luc_RNA")
samples_sirt6 <- paste0(donors, "_SIRT6_RNA")
all_samples <- c(samples_old, samples_sirt6)

# Read the gene-level TPM column from one RSEM result file.
read_sample_tpm <- function(sample, input_dir = base_dir) {
    file <- file.path(input_dir, sample, paste0(sample, ".genes.results"))
    if (!file.exists(file)) stop("Missing input file: ", file)

    result <- read.table(file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
    required_columns <- c("gene_id", "TPM")
    if (!all(required_columns %in% colnames(result))) {
        stop("Required columns are missing from: ", file)
    }

    output <- result[, required_columns]
    colnames(output)[2] <- sample
    cat("Read:", sample, "- genes:", nrow(output),
        "- TPM range:", round(min(output[[2]], na.rm = TRUE), 2), "to",
        round(max(output[[2]], na.rm = TRUE), 2), "\n")
    output
}

# Merge TPM values from all samples by gene ID.
merge_sample_tpm <- function(samples, input_dir = base_dir) {
    sample_tables <- lapply(samples, read_sample_tpm, input_dir = input_dir)
    matrix <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
                     sample_tables)
    matrix[is.na(matrix)] <- 0
    matrix
}

# Generate one annotated gene-level TPM table for all downstream analyses.
tpm_df <- merge_sample_tpm(all_samples)
stopifnot(identical(colnames(tpm_df)[-1], all_samples))

tpm_file <- file.path(output_dir, "merge.TPM.secondBatch.txt")
write.table(tpm_df, tpm_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("TPM matrix saved to:", tpm_file, "\n")
cat("Matrix dimensions:", paste(dim(tpm_df), collapse = " x "), "\n")

tpm_df <- tpm_df %>%
    mutate(ensembl = str_extract(gene_id, "^ENSG\\d+"),
           gene_symbol = mapIds(org.Hs.eg.db, keys = ensembl,
                                keytype = "ENSEMBL", column = "SYMBOL",
                                multiVals = "first"))

# Collapse multiple Ensembl IDs assigned to the same gene symbol.
expr_by_symbol <- tpm_df %>%
    filter(!is.na(gene_symbol)) %>%
    group_by(gene_symbol) %>%
    summarise(across(all_of(all_samples), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")

expr_mat <- expr_by_symbol %>%
    column_to_rownames("gene_symbol") %>%
    select(all_of(all_samples)) %>%
    as.matrix()

# Retain genes expressed above 1 TPM in every sample of at least one group.
keep_genes <- rowSums(expr_mat[, samples_old, drop = FALSE] > 1) == length(samples_old) |
              rowSums(expr_mat[, samples_sirt6, drop = FALSE] > 1) == length(samples_sirt6)
expr_mat <- expr_mat[keep_genes, , drop = FALSE]
log2_mat <- log2(expr_mat + 1)
cat("Genes retained after filtering:", nrow(log2_mat), "\n")

# Fit a paired limma model using donor as a blocking factor.
donor <- factor(rep(donors, 2))
group <- factor(c(rep("Old", length(samples_old)),
                  rep("SIRT6", length(samples_sirt6))),
                levels = c("Old", "SIRT6"))
design <- model.matrix(~ donor + group)
fit <- eBayes(lmFit(log2_mat, design))

de_results <- topTable(fit, coef = "groupSIRT6", number = Inf,
                       sort.by = "P", adjust.method = "BH") %>%
    rownames_to_column("gene_symbol")
significant_results <- de_results %>%
    filter(adj.P.Val < 0.05) %>%
    mutate(direction = ifelse(logFC > 0, "Up in SIRT6", "Down in SIRT6"))

write.csv(de_results, file.path(output_dir, "DEG_SIRT6_vs_Old_secondBatch.csv"),row.names = FALSE)
write.table(significant_results,file.path(output_dir, "DEG_SIRT6_vs_Old_significant.txt"),sep = "\t", quote = FALSE, row.names = FALSE)

cat("Significant genes (FDR < 0.05):", nrow(significant_results), "\n")
cat("Upregulated in SIRT6:", sum(significant_results$logFC > 0), "\n")
cat("Downregulated in SIRT6:", sum(significant_results$logFC < 0), "\n")

# Generate the differential-expression volcano plot.
plot_volcano <- function(results) {
    plot_data <- results %>%
        mutate(Significance = case_when(
            adj.P.Val < 0.05 & logFC > 0 ~ "Up in SIRT6",
            adj.P.Val < 0.05 & logFC < 0 ~ "Down in SIRT6",
            TRUE ~ "NS"))

    n_up <- sum(plot_data$Significance == "Up in SIRT6")
    n_down <- sum(plot_data$Significance == "Down in SIRT6")
    y_max <- max(-log10(plot_data$adj.P.Val), na.rm = TRUE)

    ggplot(plot_data, aes(logFC, -log10(adj.P.Val), color = Significance)) +
        geom_point(alpha = 0.6, size = 1.2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
        annotate("text", x = max(plot_data$logFC) * 0.9, y = y_max * 0.97,
                 label = paste0("Up: ", n_up), size = 5, hjust = 1) +
        annotate("text", x = min(plot_data$logFC) * 0.9, y = y_max * 0.97,
                 label = paste0("Down: ", n_down), size = 5, hjust = 0) +
        scale_color_manual(values = c("Up in SIRT6" = "#E74C3C",
                                      "Down in SIRT6" = "#3DBBD4",
                                      "NS" = "grey70")) +
        labs(title = "SIRT6 vs Old: limma differential expression (Batch 2)",
             x = expression(log[2]~"fold change (SIRT6/Old)"),
             y = expression(-log[10]~"FDR"), color = NULL) +
        theme_classic(base_size = 12) +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              plot.title = element_text(size = 15))
}

volcano_plot <- plot_volcano(de_results)
ggsave(file.path(output_dir, "volcano_SIRT6_vs_Old_secondBatch.pdf"),plot = volcano_plot, width = 6, height = 4)

# Generate paired expression plots for selected genes.
plot_selected_genes <- function(data, genes, old_samples, sirt6_samples) {
    selected <- data %>%
        filter(gene_symbol %in% genes) %>%
        select(gene_symbol, all_of(c(old_samples, sirt6_samples))) %>%
        pivot_longer(-gene_symbol, names_to = "sample", values_to = "TPM") %>%
        mutate(Group = ifelse(grepl("_Luc_RNA$", sample), "Old", "SIRT6"),
               Group = factor(Group, levels = c("Old", "SIRT6")),
               donor = str_remove(sample, "_Luc_RNA$|_SIRT6_RNA$"),
               gene_symbol = factor(gene_symbol, levels = genes))

    statistics <- compare_means(TPM ~ Group, data = selected,
                                group.by = "gene_symbol", method = "wilcox.test",
                                paired = TRUE) %>%
        mutate(p.adj.BH = p.adjust(p, method = "BH"),
               gene_symbol = factor(gene_symbol, levels = genes),
               p.label = paste0("p.adj = ", formatC(p.adj.BH, format = "e", digits = 2))) %>%
        left_join(selected %>% group_by(gene_symbol) %>%
                      summarise(y.position = max(TPM, na.rm = TRUE) * 1.1,
                                .groups = "drop"), by = "gene_symbol")

    ggplot(selected, aes(Group, TPM, color = Group)) +
        geom_boxplot(outlier.shape = NA, width = 0.3, fill = NA) +
        geom_line(aes(group = donor), color = "grey75", linewidth = 0.4) +
        geom_point(size = 2, alpha = 0.7) +
        geom_text(data = statistics, aes(x = 1.5, y = y.position, label = p.label),
                  inherit.aes = FALSE, size = 4) +
        scale_color_manual(values = c("Old" = "#E74C3C", "SIRT6" = "#F5A623")) +
        facet_wrap(~ gene_symbol, scales = "free_y", nrow = 1) +
        labs(y = "TPM", x = NULL, color = "Group", title = "Batch 2: Old vs SIRT6") +
        theme_classic() +
        theme(legend.position = "top",
              axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              strip.text = element_text(size = 15, face = "italic"))
}

genes_to_plot <- c("KAP1", "SUV39H1", "SUV39H2", "SETDB1", "DNMT1", "DNMT3B")
plot_data <- expr_by_symbol %>% mutate(gene_symbol = recode(gene_symbol, "TRIM28" = "KAP1"))
expression_plot <- plot_selected_genes(plot_data, genes_to_plot,samples_old, samples_sirt6)
ggsave(file.path(output_dir, "boxplot_batch2_Old_vs_SIRT6.pdf"), plot = expression_plot, width = 15, height = 4)
