# libs
library(DESeq2)
library(pheatmap)

# file paths
vsd_file <- "results/deseq2_exploration_3patients/vsd_3patients.rds"
up_file <- "results/deseq2_results_summary_3patients/upregulated_in_tumor_3patients.csv"
down_file <- "results/deseq2_results_summary_3patients/downregulated_in_tumor_3patients.csv"
output_dir <- "results/deseq2_top_heatmap_3patients"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# number of genes to keep
top_n_up <- 20
top_n_down <- 20

# necessary files load in
vsd <- readRDS(vsd_file)
up_genes <- read.csv(up_file, stringsAsFactors = FALSE)
down_genes <- read.csv(down_file, stringsAsFactors = FALSE)

# keeps only genes with symbols
up_genes <- up_genes[!is.na(up_genes$gene_symbol) & up_genes$gene_symbol != "", , drop = FALSE]
down_genes <- down_genes[!is.na(down_genes$gene_symbol) & down_genes$gene_symbol != "", , drop = FALSE]

# keeps top genes after symbol filtering
top_up <- head(up_genes, top_n_up)
top_down <- head(down_genes, top_n_down)

top_genes <- rbind(top_up, top_down)

# extracts matrix to be used for heatmap
vst_mat <- assay(vsd)

# filters gene that match vst matrix
selected_ids <- top_genes$ensembl_gene_id
selected_ids <- selected_ids[selected_ids %in% rownames(vst_mat)]

heatmap_mat <- vst_mat[selected_ids, , drop = FALSE]

# matches top gene annotation back to the heatmap amtrix order
top_genes <- top_genes[match(rownames(heatmap_mat), top_genes$ensembl_gene_id), ,drop = FALSE]

# row labels
top_genes$plot_label <- top_genes$gene_symbol
top_genes$plot_label <- make.unique(top_genes$plot_label)

rownames(heatmap_mat) <- top_genes$plot_label

# scales each gene across samples
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# removes rows of failed scaling
keep_rows <- complete.cases(heatmap_mat_scaled)
heatmap_mat_scaled <- heatmap_mat_scaled[keep_rows,, drop = FALSE]
top_genes <- top_genes[keep_rows, , drop = FALSE]

# sample annotation from vst
sample_info <- as.data.frame(colData(vsd))
sample_info$sample_name <- rownames(sample_info)

sample_info$condition <- factor(sample_info$condition, levels = c("normal", "tumor"))
sample_info$patient <- factor(
  as.character(sample_info$patient),
  levels = sort(unique(as.character(sample_info$patient)))
)

sample_info <- sample_info[order(sample_info$patient, sample_info$condition), , drop = FALSE]

# reorders matrix collumns to match annotation
heatmap_mat_scaled <- heatmap_mat_scaled[,rownames(sample_info), drop = FALSE]

annotation_col <- sample_info[, c("condition", "patient"), drop = FALSE]

# data saving
write.csv(
  top_genes,
  file = file.path(output_dir, "top_genes_used_for_heatmap_3patients.csv"),
  row.names = FALSE
)

write.csv(
  heatmap_mat_scaled,
  file = file.path(output_dir, "heatmap_matrix_scaled_3patients.csv"),
  row.names = FALSE
)

write.csv(
  sample_info,
  file = file.path(output_dir,"heatmap_sample_info_3patients.csv"),
  row.names = TRUE
)

#heatmap
heatmap_mat_plot <- heatmap_mat_scaled

# shorten column labels
colnames(heatmap_mat_plot) <- c("2_N", "2_T", "3_N", "3_T", "5_N", "5_T")

# keeps only condition annotation
annotation_col_plot <- data.frame(
  condition = sample_info$condition
)

rownames(annotation_col_plot) <- colnames(heatmap_mat_plot)

# annotation colors
annotation_colors <- list(
  condition = c(
    normal = "#9cf4a5",
    tumor = "#D97CE8"
  )
)

# heatmap colors
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# save 
pheatmap(
  heatmap_mat_plot,
  color = heatmap_colors,
  annotation_col = annotation_col_plot,
  annotation_colors = annotation_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = c(2, 4),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 11,
  angle_col = 0,
  border_color = NA,
  main = "Top differentially expressed genes: tumor vs normal",
  filename = file.path(output_dir, "heatmap_top_de_genes_3patients.png"),
  width = 7,
  height = 10
)

paste("heatmap created in: ",output_dir)