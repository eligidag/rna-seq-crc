# libraries
library(tximport)

# files
analysis_subset_file <- "data/metadata/gse50760_analysis_subset_3patients.tsv"
tx2gene_file <- "reference/gencode_v49_salmon/tx2gene.tsv"
salmon_dir <- "results/salmon_quant_3patients/"
output_dir <- "results/tximport_3patients"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_table <- read.delim(
  file = analysis_subset_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# condition label to be reused later
sample_table$condition <- ifelse(
  sample_table$sample_type == "primary colorectal cancer",
  "tumor",
  "normal"
)

# rebuilds sample names to match the Salmon output folders
sample_table$sample_name <- paste0(
  "AMC_",
  sample_table$patient_id,
  "_",
  sample_table$condition
)

sample_table$quant_file <- file.path(
  salmon_dir,
  sample_table$sample_name,
  "quant.sf"
)

# stops early if any of the expected Salomn results are missing
missing_quant_files <- sample_table$quant_file[!file.exists(sample_table$quant_file)]

if(length(missing_quant_files) > 0) {
  stop(
    paste(
      "stopped, expected quant.sf file/s are missing:",
      paste(missing_quant_files, collapse =", ")
    )
  )
}

# reads the tx2gene mapping table
tx2gene <- read.delim(
  file = tx2gene_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# builds the named vector of Salmon quant files
files <-sample_table$quant_file
names(files) <- sample_table$sample_name

# imports Salmon quantification and summarizes to gene level
txi <- tximport(
  files = files,
  type = "salmon",
  tx2gene = tx2gene
)

# rds
saveRDS(
  object = txi,
  file = file.path(output_dir, "txi_3patients.rds")
)

# gene counts
write.csv(
  x = as.data.frame(txi$counts),
  file = file.path(output_dir, "gene_counts_3patients.csv")
)

# gene length
write.csv(
  x = as.data.frame(txi$length),
  file = file.path(output_dir, "gene_length_3patients.csv")
)

write.csv(
  x = sample_table,
  file = file.path(output_dir, "sample_table_3patients.csv"),
  row.names = FALSE
)

write.csv(
  x = as.data.frame(txi$abundance),
  file = file.path(output_dir, "gene_abundance_tpm_3patients.csv")
)

# summary
cat("tximport completed successfully\n")
cat("Number of samples imported:", ncol(txi$counts), "\n")
cat("Number of genes summarised:", nrow(txi$counts), "\n")
cat("Results written to:", output_dir, "\n")
