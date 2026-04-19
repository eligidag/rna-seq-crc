#libs
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

# file paths
results_rds_file <- "results/deseq2_3patients/deseq2_results_3patients.rds"
output_dir <- "results/deseq2_results_summary_3patients/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1

# loads DESeq2 results
res <- readRDS(results_rds_file)

# converts results object into workable df
res_df <- as.data.frame(res)

# extracts Ensembl genes ids into a columm
res_df$ensembl_gene_id  <- rownames(res_df)

# removes version numbers from Ensembl ids
res_df$ensembl_gene_id_clean <- sub("\\..*", "", res_df$ensembl_gene_id)

# maps Ensembl ids to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ensembl_gene_id_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# adds symbols to the df
res_df$gene_symbol <- unname(gene_symbols)
res_df$gene_symbol[is.na(res_df$gene_symbol)] <- ""

# sorts df based on relevance
res_df <- res_df[, c(
  "ensembl_gene_id",
  "ensembl_gene_id_clean",
  "gene_symbol",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

# removes rows with missing adjusted p-values
res_df_tested <- res_df[!is.na(res_df$padj),]

# classifies genes by significance and direction
res_df_tested$regulation <- "not_significant"

res_df_tested$regulation[
  res_df_tested$padj < padj_cutoff & res_df_tested$log2FoldChange >= log2fc_cutoff
] <- "up_in_tumor"

res_df_tested$regulation[
  res_df_tested$padj < padj_cutoff & res_df_tested$log2FoldChange <= -log2fc_cutoff
] <- "down_in_tumor"

# splits significant genes into individual tables
sig_genes <- res_df_tested[res_df_tested$padj < padj_cutoff,]

up_genes <- res_df_tested[res_df_tested$regulation == "up_in_tumor",]
down_genes <- res_df_tested[res_df_tested$regulation == "down_in_tumor",]

# sorts tables of interest by adjusted p-value
sig_genes <- sig_genes[order(sig_genes$padj),]

up_genes <- up_genes[order(up_genes$padj),]
down_genes <- down_genes[order(down_genes$padj),]

# top 20 up/down regulated genes
top_20_up <- head(up_genes,20) 
top_20_down <- head(down_genes,20)

# saves results for future processing
write.csv(res_df,
          file = file.path(output_dir, "all_results_cleaned_3patients.csv"),
          row.names = FALSE)

write.csv(sig_genes,
          file = file.path(output_dir, "significant_genes_3patients.csv"),
          row.names = FALSE)

write.csv(up_genes,
          file = file.path(output_dir, "upregulated_in_tumor_3patients.csv"),
          row.names = FALSE)

write.csv(down_genes,
          file = file.path(output_dir, "downregulated_in_tumor_3patients.csv"),
          row.names = FALSE)

write.csv(top_20_up,
          file = file.path(output_dir, "top_20_upregulated_in_tumor_3patients.csv"),
          row.names = FALSE)

write.csv(top_20_down,
          file = file.path(output_dir, "top_20_downregulated_in_tumor_3patients.csv"),
          row.names = FALSE)

# numbers for summary
n_all <- nrow(res_df)
n_tested <- nrow(res_df_tested)
n_sig <- nrow(sig_genes)
n_up <- nrow(up_genes)
n_down <- nrow(down_genes)

#summary
summary_lines <- c(
  "DESeq2 results summary",
  paste("Total genes in results table:", n_all),
  paste("Genes with adjusted p-values:", n_tested),
  paste("Significant genes (padj < 0.05):", n_sig),
  paste("Upregulated genes in tumor (padj < 0.05 and log2FC >= 1):",n_up),
  paste("Downregulated genes in tumor (padj < 0.05 and log2FC <= -1):", n_down)
)

summary_lines

writeLines(summary_lines,
           con = file.path(output_dir, "summary_stats_3patients.txt"))

# prep for volcano plot
res_df_tested$neg_log10_padj <- -log10(res_df_tested$padj)

# volcano plot
volcano_plot <- ggplot(
  res_df_tested,
  aes(x = log2FoldChange, y= neg_log10_padj, color = regulation)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c(
      "not_significant" = "grey70",
      "up_in_tumor" = "blue",
      "down_in_tumor" = "red"
    )
  ) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  labs(
    title = "Volcano plot: tumor vs normal",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value",
    color = "Gene group"
  ) + 
  theme_bw()

ggsave(
  filename = file.path(output_dir, "volcano_plot_3patients.png"),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

paste("Volcano plot saved in:",output_dir)


# picks top 5 upregulated and downregulated genes for labeling
top_5_up_labels <- head(up_genes, 5)
top_5_down_labels <- head(down_genes, 5)

label_df <- rbind(top_5_up_labels, top_5_down_labels)

# y axis values needed for labeling
label_df$neg_log10_padj <- -log10(label_df$padj)

# uses gene symbol when available
label_df$plot_label <- label_df$gene_symbol
label_df$plot_label[label_df$plot_label == ""] <- label_df$ensembl_gene_id_clean[label_df$plot_label == ""]

#plot
volcano_plot_labeled <- ggplot(
  res_df_tested,
  aes(x = log2FoldChange, y = neg_log10_padj, color = regulation)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c(
      "not_significant" = "grey70",
      "up_in_tumor" = "blue",
      "down_in_tumor" = "red"
    )
  ) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = log2FoldChange, y = neg_log10_padj, label = plot_label),
    size = 3.5,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  labs(
    title = "Volcano plot: tumor vs normal (top 5 up/down)",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value",
    color = "Gene group"
  ) +
  theme_bw()


ggsave(
  filename = file.path(output_dir, "volcano_plot_3patients_top5_up_down.png"),
  plot = volcano_plot_labeled,
  width = 8,
  height = 6,
  dpi = 300
)

paste("Volcano plot with top 5 up/down regulated genes labels saved in:",output_dir)