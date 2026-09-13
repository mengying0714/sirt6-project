#!/bin/bash

module load bluehive/2.5
module load r/4.2.1/b1

BASE_DIR="/RNA_seq/STAR_New"
OUTPUT_DIR="${BASE_DIR}/count1"
mkdir -p "${OUTPUT_DIR}"

declare -A sample_group
sample_group["61F_Luc_RNA"]="Old"
sample_group["61F_SIRT6_RNA"]="SIRT6"
sample_group["70F_Luc_RNA"]="Old"
sample_group["70F_SIRT6_RNA"]="SIRT6"
sample_group["71F_Luc_RNA"]="Old"
sample_group["71F_SIRT6_RNA"]="SIRT6"
sample_group["73M_Luc_RNA"]="Old"
sample_group["73M_SIRT6_RNA"]="SIRT6"
sample_group["74F_Luc_RNA"]="Old"
sample_group["74F_SIRT6_RNA"]="SIRT6"
sample_group["74M_Luc_RNA"]="Old"
sample_group["74M_SIRT6_RNA"]="SIRT6"
sample_group["79M_Luc_RNA"]="Old"
sample_group["79M_SIRT6_RNA"]="SIRT6"
sample_group["85F_Luc_RNA"]="Old"
sample_group["85F_SIRT6_RNA"]="SIRT6"
sample_group["86F_Luc_RNA"]="Old"
sample_group["86F_SIRT6_RNA"]="SIRT6"
sample_group["88M_Luc_RNA"]="Old"
sample_group["88M_SIRT6_RNA"]="SIRT6"

for sample in "${!sample_group[@]}"; do
    results_file="${BASE_DIR}/${sample}/${sample}.genes.results"
    if [[ -f "${results_file}" ]]; then
        echo "Processing ${sample} (${sample_group[$sample]})..."
        awk 'NR > 1 {print $1"\t"$5}' "${results_file}" > \
            "${OUTPUT_DIR}/${sample}.counts.txt"
    else
        echo "Warning: ${results_file} not found!"
    fi
done

echo "sample,condition,donor" > "${OUTPUT_DIR}/sample_info_batch2.csv"
for sample in "${!sample_group[@]}"; do
    donor=$(echo "${sample}" | sed 's/_Luc_RNA//;s/_SIRT6_RNA//')
    echo "${sample},${sample_group[$sample]},${donor}"
done >> "${OUTPUT_DIR}/sample_info_batch2.csv"

echo "Extraction complete. Starting data merging in R..."

cat > "${OUTPUT_DIR}/combine_counts_batch2.R" << 'EOF'
library(dplyr)

samples_old <- c(
    "61F_Luc_RNA", "70F_Luc_RNA", "71F_Luc_RNA", "73M_Luc_RNA",
    "74F_Luc_RNA", "74M_Luc_RNA", "79M_Luc_RNA", "85F_Luc_RNA",
    "86F_Luc_RNA", "88M_Luc_RNA"
)

samples_sirt6 <- c(
    "61F_SIRT6_RNA", "70F_SIRT6_RNA", "71F_SIRT6_RNA", "73M_SIRT6_RNA",
    "74F_SIRT6_RNA", "74M_SIRT6_RNA", "79M_SIRT6_RNA", "85F_SIRT6_RNA",
    "86F_SIRT6_RNA", "88M_SIRT6_RNA"
)

all_samples <- c(samples_old, samples_sirt6)

count_list <- lapply(all_samples, function(s) {
    f <- paste0(s, ".counts.txt")

    if (!file.exists(f)) {
        cat("Missing:", f, "\n")
        return(NULL)
    }

    df <- read.table(
        f,
        header = FALSE,
        col.names = c("gene_id", s),
        stringsAsFactors = FALSE
    )

    cat("Read:", s, "- genes:", nrow(df), "\n")
    df
})

# Stop if one or more count files are missing
if (any(vapply(count_list, is.null, logical(1)))) {
    stop("One or more count files are missing.")
}

count_matrix <- Reduce(
    function(x, y) merge(x, y, by = "gene_id", all = TRUE),
    count_list
)

count_matrix[is.na(count_matrix)] <- 0
rownames(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix[, -1, drop = FALSE]

cat("\nColumn-name verification:\n")
cat("Old samples:", colnames(count_matrix)[1:10], "\n")
cat("SIRT6 samples:", colnames(count_matrix)[11:20], "\n")

# Verify that the sample columns are correctly ordered
cat("\nConfirm that all Old columns contain 'Luc':\n")
print(all(grepl("Luc", colnames(count_matrix)[1:10])))

cat("Confirm that all SIRT6 columns contain 'SIRT6':\n")
print(all(grepl("SIRT6", colnames(count_matrix)[11:20])))

write.table(
    count_matrix,
    file = "secondBatch_count_matrix_correct.txt",
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
)

sample_info <- data.frame(
    sample = all_samples,
    condition = c(rep("Old", 10), rep("SIRT6", 10)),
    donor = gsub("_Luc_RNA|_SIRT6_RNA", "", all_samples),
    row.names = all_samples
)

write.csv(
    sample_info,
    "sample_info_batch2_correct.csv",
    row.names = FALSE
)

cat("\nSample information:\n")
print(sample_info)
cat("\nComplete! Matrix dimensions:", dim(count_matrix), "\n")
EOF

cd "${OUTPUT_DIR}" || exit 1
Rscript combine_counts_batch2.R
echo "Complete!"