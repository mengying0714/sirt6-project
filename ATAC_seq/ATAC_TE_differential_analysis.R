###############################################################################
# ATAC-seq RepeatMasker RPKM calculation and paired differential analysis
###############################################################################

library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(Rmisc)

# Set directories through environment variables.
data_dir <- Sys.getenv("ATAC_DATA_DIR", unset = ".")
output_dir <- Sys.getenv("ATAC_OUTPUT_DIR", unset = ".")
annotation_dir <- Sys.getenv("ATAC_ANNOTATION_DIR", unset = ".")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Part 1: Read RepeatMasker counts and calculate RPKM
###############################################################################

S61F_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.61F_Luc.Reads.Num.txt"),header = TRUE, sep = "")
S71F_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.71F_Luc.Reads.Num.txt"),header = FALSE, sep = "")
S73M_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.73M_Luc.Reads.Num.txt"),header = FALSE, sep = "")
S74F_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.74F_Luc.Reads.Num.txt"),header = FALSE, sep = "")
S85F_Luc <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.85F_Luc.Reads.Num.txt"),header = FALSE, sep = "")

S61F_SIRT6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.61F_SIRT6.Reads.Num.txt"),header = FALSE, sep = "")
S71F_SIRT6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.71F_SIRT6.Reads.Num.txt"),header = FALSE, sep = "")
S73M_SIRT6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.73M_SIRT6.Reads.Num.txt"),header = FALSE, sep = "")
S74F_SIRT6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.74F_SIRT6.Reads.Num.txt"),header = FALSE, sep = "")
S85F_SIRT6 <- read.table(file.path(data_dir, "RepeatMasker.hg38.Family.repName.85F_SIRT6.Reads.Num.txt"),header = FALSE, sep = "")

# Assign consistent column names.
colnames(S61F_Luc) <- c("repName", "ReadNum")
colnames(S71F_Luc) <- c("repName", "ReadNum")
colnames(S73M_Luc) <- c("repName", "ReadNum")
colnames(S74F_Luc) <- c("repName", "ReadNum")
colnames(S85F_Luc) <- c("repName", "ReadNum")

colnames(S61F_SIRT6) <- c("repName", "ReadNum")
colnames(S71F_SIRT6) <- c("repName", "ReadNum")
colnames(S73M_SIRT6) <- c("repName", "ReadNum")
colnames(S74F_SIRT6) <- c("repName", "ReadNum")
colnames(S85F_SIRT6) <- c("repName", "ReadNum")

# Replace missing read counts with zero.
S61F_Luc[is.na(S61F_Luc)] <- 0
S71F_Luc[is.na(S71F_Luc)] <- 0
S73M_Luc[is.na(S73M_Luc)] <- 0
S74F_Luc[is.na(S74F_Luc)] <- 0
S85F_Luc[is.na(S85F_Luc)] <- 0

S61F_SIRT6[is.na(S61F_SIRT6)] <- 0
S71F_SIRT6[is.na(S71F_SIRT6)] <- 0
S73M_SIRT6[is.na(S73M_SIRT6)] <- 0
S74F_SIRT6[is.na(S74F_SIRT6)] <- 0
S85F_SIRT6[is.na(S85F_SIRT6)] <- 0

# Read the total genomic length of each RepeatMasker repName.
RepName.Length <- read.table(file.path(annotation_dir,"RepeatMasker.hg38.Family.repName.merged.sorted.Total.Length.txt"),header = TRUE,sep = "\t")

# Add repeat-family length information to each sample.
S61F_Luc.Length <- right_join(S61F_Luc, RepName.Length, by = "repName")
S71F_Luc.Length <- right_join(S71F_Luc, RepName.Length, by = "repName")
S73M_Luc.Length <- right_join(S73M_Luc, RepName.Length, by = "repName")
S74F_Luc.Length <- right_join(S74F_Luc, RepName.Length, by = "repName")
S85F_Luc.Length <- right_join(S85F_Luc, RepName.Length, by = "repName")

S61F_SIRT6.Length <- right_join(S61F_SIRT6, RepName.Length, by = "repName")
S71F_SIRT6.Length <- right_join(S71F_SIRT6, RepName.Length, by = "repName")
S73M_SIRT6.Length <- right_join(S73M_SIRT6, RepName.Length, by = "repName")
S74F_SIRT6.Length <- right_join(S74F_SIRT6, RepName.Length, by = "repName")
S85F_SIRT6.Length <- right_join(S85F_SIRT6, RepName.Length, by = "repName")

# Add the total number of duplicate-removed reads for each sample.
S61F_Luc.Length$TotalRead <- 138120413
S71F_Luc.Length$TotalRead <- 112703209
S73M_Luc.Length$TotalRead <- 104366319
S74F_Luc.Length$TotalRead <- 97842102
S85F_Luc.Length$TotalRead <- 94798565

S61F_SIRT6.Length$TotalRead <- 80824640
S71F_SIRT6.Length$TotalRead <- 95911784
S73M_SIRT6.Length$TotalRead <- 104750833
S74F_SIRT6.Length$TotalRead <- 131907666
S85F_SIRT6.Length$TotalRead <- 62596670

# Convert read counts to numeric values.
S61F_Luc.Length$ReadNum <- as.numeric(S61F_Luc.Length$ReadNum)
S71F_Luc.Length$ReadNum <- as.numeric(S71F_Luc.Length$ReadNum)
S73M_Luc.Length$ReadNum <- as.numeric(S73M_Luc.Length$ReadNum)
S74F_Luc.Length$ReadNum <- as.numeric(S74F_Luc.Length$ReadNum)
S85F_Luc.Length$ReadNum <- as.numeric(S85F_Luc.Length$ReadNum)

S61F_SIRT6.Length$ReadNum <- as.numeric(S61F_SIRT6.Length$ReadNum)
S71F_SIRT6.Length$ReadNum <- as.numeric(S71F_SIRT6.Length$ReadNum)
S73M_SIRT6.Length$ReadNum <- as.numeric(S73M_SIRT6.Length$ReadNum)
S74F_SIRT6.Length$ReadNum <- as.numeric(S74F_SIRT6.Length$ReadNum)
S85F_SIRT6.Length$ReadNum <- as.numeric(S85F_SIRT6.Length$ReadNum)

# Calculate RPKM for each RepeatMasker repName.
S61F_Luc.Length$S61F_Luc.RPKM <- S61F_Luc.Length$ReadNum /S61F_Luc.Length$TotalRead / S61F_Luc.Length$TotalLength * 1e9

S71F_Luc.Length$S71F_Luc.RPKM <- S71F_Luc.Length$ReadNum /S71F_Luc.Length$TotalRead / S71F_Luc.Length$TotalLength * 1e9

S73M_Luc.Length$S73M_Luc.RPKM <- S73M_Luc.Length$ReadNum /S73M_Luc.Length$TotalRead / S73M_Luc.Length$TotalLength * 1e9

S74F_Luc.Length$S74F_Luc.RPKM <- S74F_Luc.Length$ReadNum /S74F_Luc.Length$TotalRead / S74F_Luc.Length$TotalLength * 1e9

S85F_Luc.Length$S85F_Luc.RPKM <- S85F_Luc.Length$ReadNum /S85F_Luc.Length$TotalRead / S85F_Luc.Length$TotalLength * 1e9

S61F_SIRT6.Length$S61F_SIRT6.RPKM <- S61F_SIRT6.Length$ReadNum /S61F_SIRT6.Length$TotalRead / S61F_SIRT6.Length$TotalLength * 1e9

S71F_SIRT6.Length$S71F_SIRT6.RPKM <- S71F_SIRT6.Length$ReadNum /S71F_SIRT6.Length$TotalRead / S71F_SIRT6.Length$TotalLength * 1e9

S73M_SIRT6.Length$S73M_SIRT6.RPKM <- S73M_SIRT6.Length$ReadNum /S73M_SIRT6.Length$TotalRead / S73M_SIRT6.Length$TotalLength * 1e9

S74F_SIRT6.Length$S74F_SIRT6.RPKM <- S74F_SIRT6.Length$ReadNum /S74F_SIRT6.Length$TotalRead / S74F_SIRT6.Length$TotalLength * 1e9

S85F_SIRT6.Length$S85F_SIRT6.RPKM <- S85F_SIRT6.Length$ReadNum /S85F_SIRT6.Length$TotalRead / S85F_SIRT6.Length$TotalLength * 1e9

# Combine the RPKM values from all samples.
ATAC.RPKM <- left_join(
  S61F_Luc.Length[, c("repName", "S61F_Luc.RPKM")],
  S71F_Luc.Length[, c("repName", "S71F_Luc.RPKM")],
  by = "repName"
)
ATAC.RPKM <- left_join(
  ATAC.RPKM,
  S73M_Luc.Length[, c("repName", "S73M_Luc.RPKM")],
  by = "repName"
)
ATAC.RPKM <- left_join(
  ATAC.RPKM,
  S74F_Luc.Length[, c("repName", "S74F_Luc.RPKM")],
  by = "repName"
)
ATAC.RPKM <- left_join(
  ATAC.RPKM,
  S85F_Luc.Length[, c("repName", "S85F_Luc.RPKM")],
  by = "repName"
)
ATAC.RPKM <- left_join(
  ATAC.RPKM,
  S61F_SIRT6.Length[, c("repName", "S61F_SIRT6.RPKM")],
  by = "repName"
)
ATAC.RPKM <- left_join(
  ATAC.RPKM,
  S71F_SIRT6.Length[, c("repName", "S71F_SIRT6.RPKM")],
  by = "repName"
)
ATAC.RPKM <- left_join(
  ATAC.RPKM,
  S73M_SIRT6.Length[, c("repName", "S73M_SIRT6.RPKM")],
  by = "repName"
)
ATAC.RPKM <- left_join(
  ATAC.RPKM,
  S74F_SIRT6.Length[, c("repName", "S74F_SIRT6.RPKM")],
  by = "repName"
)
ATAC.RPKM <- left_join(
  ATAC.RPKM,
  S85F_SIRT6.Length[, c("repName", "S85F_SIRT6.RPKM")],
  by = "repName"
)

# Combine the raw read counts from all samples.
ATAC.Count <- left_join(
  S61F_Luc.Length[, c("repName", "ReadNum")],
  S71F_Luc.Length[, c("repName", "ReadNum")],
  by = "repName"
)
ATAC.Count <- left_join(
  ATAC.Count,
  S73M_Luc.Length[, c("repName", "ReadNum")],
  by = "repName"
)
ATAC.Count <- left_join(
  ATAC.Count,
  S74F_Luc.Length[, c("repName", "ReadNum")],
  by = "repName"
)
ATAC.Count <- left_join(
  ATAC.Count,
  S85F_Luc.Length[, c("repName", "ReadNum")],
  by = "repName"
)
ATAC.Count <- left_join(
  ATAC.Count,
  S61F_SIRT6.Length[, c("repName", "ReadNum")],
  by = "repName"
)
ATAC.Count <- left_join(
  ATAC.Count,
  S71F_SIRT6.Length[, c("repName", "ReadNum")],
  by = "repName"
)
ATAC.Count <- left_join(
  ATAC.Count,
  S73M_SIRT6.Length[, c("repName", "ReadNum")],
  by = "repName"
)
ATAC.Count <- left_join(
  ATAC.Count,
  S74F_SIRT6.Length[, c("repName", "ReadNum")],
  by = "repName"
)
ATAC.Count <- left_join(
  ATAC.Count,
  S85F_SIRT6.Length[, c("repName", "ReadNum")],
  by = "repName"
)

colnames(ATAC.Count) <- c(
  "repName",
  "S61F_Luc", "S71F_Luc", "S73M_Luc", "S74F_Luc", "S85F_Luc",
  "S61F_SIRT6", "S71F_SIRT6", "S73M_SIRT6",
  "S74F_SIRT6", "S85F_SIRT6"
)

colnames(ATAC.RPKM) <- c(
  "repName",
  "S61F_Luc", "S71F_Luc", "S73M_Luc", "S74F_Luc", "S85F_Luc",
  "S61F_SIRT6", "S71F_SIRT6", "S73M_SIRT6",
  "S74F_SIRT6", "S85F_SIRT6"
)

save(
  ATAC.Count,
  ATAC.RPKM,
  file = file.path(output_dir, "ATAC.Res.RData")
)

###############################################################################
# Part 2: Paired differential ATAC-seq analysis using edgeR
###############################################################################

ATAC.Count[is.na(ATAC.Count)] <- 0
ATAC.Count <- na.omit(ATAC.Count)

rownames(ATAC.Count) <- ATAC.Count$repName
rownames(ATAC.RPKM) <- ATAC.RPKM$repName

samples <- read.table(
  file.path(data_dir, "ATAC.HDF.samples.txt"),
  header = TRUE,
  sep = "\t"
)

TE.infor <- read.table(
  file.path(annotation_dir, "TE.infor.txt"),
  header = TRUE,
  sep = "\t"
)

y <- DGEList(
  counts = ATAC.Count[, 2:11],
  group = samples$Group,
  lib.size = samples$lib.size
)

# Normalize library sizes.
y <- calcNormFactors(y)

# Generate MDS plots before dispersion estimation.
pdf(
  file.path(output_dir, "ATAC_MDS_plot.pdf"),
  width = 6,
  height = 6
)
plotMDS(y)
plotMDS(y, labels = y$samples$group)
dev.off()

# Estimate dispersion using the paired donor design.
design <- model.matrix(~ subject + Group, data = samples)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

# Test the SIRT6 versus Luc group coefficient.
group_coef <- grep("^Group", colnames(design))

if (length(group_coef) != 1) {
  stop(
    "The group coefficient could not be identified. Check samples$Group ",
    "and the design matrix."
  )
}

qlf <- glmQLFTest(fit, coef = group_coef)

# Extract all differential results.
top <- topTags(
  qlf,
  n = nrow(y$counts)
)$table

TE.infor1 <- TE.infor[rownames(top), ]
RPKM <- ATAC.RPKM[rownames(top), ]

edgeR_qlf_RPKM <- cbind(
  TE.infor1,
  top,
  RPKM
)

# Classify differential ATAC signals using the original P-value threshold.
edgeR_qlf_RPKM$ATAC.Diff <- "2.Not"

edgeR_qlf_RPKM$ATAC.Diff[
  edgeR_qlf_RPKM$PValue < 0.05 &
    edgeR_qlf_RPKM$logFC > 0
] <- "1.UP"

edgeR_qlf_RPKM$ATAC.Diff[
  edgeR_qlf_RPKM$PValue < 0.05 &
    edgeR_qlf_RPKM$logFC < 0
] <- "0.Down"

write.table(
  edgeR_qlf_RPKM,
  file.path(output_dir, "edgeR_ATAC_Paired.txt"),
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

edgeR_ATAC <- na.omit(edgeR_qlf_RPKM)

edgeR_ATAC.LINE <- subset(edgeR_ATAC, Family == "LINE")
edgeR_ATAC.SINE <- subset(edgeR_ATAC, Family == "SINE")
edgeR_ATAC.LTR <- subset(edgeR_ATAC, Family == "LTR")
edgeR_ATAC.DNA <- subset(edgeR_ATAC, Family == "DNA")

save(
  edgeR_ATAC,
  ATAC.Count,
  ATAC.RPKM,
  file = file.path(output_dir, "ATAC.Res.RData")
)

###############################################################################
# Part 3: Volcano plots of differential ATAC signals
###############################################################################

pdf(
  file = file.path(
    output_dir,
    "Volcano_Plot_Differential_ATAC_Paired.pdf"
  ),
  width = 12,
  height = 3
)

p1 <- ggplot(
  edgeR_ATAC.LINE,
  aes(x = logFC, y = -log10(PValue))
) +
  geom_point(aes(color = ATAC.Diff), size = 1.5) +
  scale_color_manual(
    values = c(
      "0.Down" = "blue",
      "1.UP" = "red",
      "2.Not" = "black"
    )
  ) +
  ggtitle("LINE") +
  xlim(-0.4, 0.4) +
  ylim(0, 3) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

p2 <- ggplot(
  edgeR_ATAC.SINE,
  aes(x = logFC, y = -log10(PValue))
) +
  geom_point(aes(color = ATAC.Diff), size = 1.5) +
  scale_color_manual(
    values = c(
      "0.Down" = "blue",
      "1.UP" = "red",
      "2.Not" = "black"
    )
  ) +
  ggtitle("SINE") +
  xlim(-0.2, 0.2) +
  ylim(0, 1.5) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

p3 <- ggplot(
  edgeR_ATAC.LTR,
  aes(x = logFC, y = -log10(PValue))
) +
  geom_point(aes(color = ATAC.Diff), size = 1.5) +
  scale_color_manual(
    values = c(
      "0.Down" = "blue",
      "1.UP" = "red",
      "2.Not" = "black"
    )
  ) +
  ggtitle("LTR") +
  xlim(-0.55, 0.55) +
  ylim(0, 4) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

p4 <- ggplot(
  edgeR_ATAC.DNA,
  aes(x = logFC, y = -log10(PValue))
) +
  geom_point(aes(color = ATAC.Diff), size = 1.5) +
  scale_color_manual(
    values = c(
      "0.Down" = "blue",
      "1.UP" = "red",
      "2.Not" = "black"
    )
  ) +
  ggtitle("DNA transposon") +
  xlim(-0.6, 0.6) +
  ylim(0, 4) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

multiplot(p1, p2, p3, p4, cols = 4)
dev.off()