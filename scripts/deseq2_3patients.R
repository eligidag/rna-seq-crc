# libs

library(DESeq2)

# file paths
txi_file <- "results/tximport_3patients/txi_3patients.rds"
metadata_file <- "data/metadata/gse50760_analysis_subset_3patients.tsv"
output_dir <- "results/deseq2_3patients"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# loads tximport output and sample metadata
txi <- readRDS(txi_file)

sample_table <- read.delim(
  file = metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# builds condition label
sample_table$condition <- ifelse(
  sample_table$sample_type == "primary colorectal cancer",
  "tumor",
  "normal"
)

# rebuilds sample names to match the tximport column names
sample_table$sample_name <- paste0(
  "AMC_",
  sample_table$patient_id,
  "_",
  sample_table$condition
)

# saves columns needed for DESeq2
coldata <- sample_table[, c("sample_name", "patient_id", "condition")]
colnames(coldata) <- c("sample", "patient", "condition")

# reorders metadata to match the count matrix
coldata <- coldata[match(colnames(txi$counts), coldata$sample),]

# stops the script if matching fails
if (any(is.na(coldata$sample))) {
  stop("matching failed between tximport counts and metadata")
}

rownames(coldata) <- coldata$sample

# converts grouping variables to factors
coldata$patient <- factor(coldata$patient)
coldata$condition <- factor(coldata$condition, levels = c("normal", "tumor"))

# builds DESeq2 object
dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = coldata,
  design = ~ patient + condition
)

# removes low count genes
dds <- dds[rowSums(counts(dds)) > 1,]

# runs DESeq2
dds <- DESeq(dds)

# tumor vs normal results
res <- results(
  dds,
  contrast = c("condition", "tumor", "normal")
)

# orders results based on adjusted p-value
res_ordered <- res[order(res$padj),]

# file saves

saveRDS(
  object = dds,
  file = file.path(output_dir, "dds_3patients.rds")
)

saveRDS(
  object = res,
  file = file.path(output_dir, "deseq2_results_3patients.rds")
)

# metadata
write.csv(
  x = coldata,
  file = file.path(output_dir, "coldata_3patients.csv"),
  row.names = FALSE
)

# normalised counts
write.csv(
  x = as.data.frame(counts(dds, normalized = TRUE)),
  file = file.path(output_dir, "normalized_counts_3patients.csv")
)

# ordered results 
write.csv(
  x = as.data.frame(res_ordered),
  file = file.path(output_dir, "deseq2_results_ordered_3patients.csv")
)

#summary
cat("DESeq2 completed successfully\n")
cat("Number of samples:", ncol(dds), "\n")
cat("Number of genes after filtering:", nrow(dds), "\n")
cat("Results written to:", output_dir, "\n")
