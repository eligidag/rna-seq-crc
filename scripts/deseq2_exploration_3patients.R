# libs
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# file paths
dds_file <- "results/deseq2_3patients/dds_3patients.rds"
output_dir <- "results/deseq2_exploration_3patients"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# DESeq2 object
dds <- readRDS(dds_file)

# makes a variance stabilised version of the count data for PCA and clustering
vsd <- vst(dds, blind = FALSE)

saveRDS(
  object = vsd,
  file = file.path(output_dir, "vsd_3patients.rds")
)

# extracts transformer expression matrix
vsd_mat <- assay(vsd)

write.csv(
  x = as.data.frame(vsd_mat),
  file = file.path(output_dir, "vst_matrix_3patients.csv")
)

# metadata table from the DESeq2 object
sample_info <- as.data.frame(colData(vsd))
sample_info$sample <- rownames(sample_info)

write.csv(
  x = sample_info,
  file = file.path(output_dir, "sample_info_3patients.csv"),
  row.names = FALSE
)

# runs PCA on the transformed data 
pca <- prcomp(t(vsd_mat), scale. = FALSE)

# df for plotting
pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)
pca_df$patient <- sample_info$patient[match(pca_df$sample, sample_info$sample)]
pca_df$condition <- sample_info$condition[match(pca_df$sample, sample_info$sample)]

# calculates percent variance
percent_var <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_var <- round(percent_var[1] * 100, 1)
pc2_var <- round(percent_var[2] * 100, 1)

write.csv(
  x = pca_df,
  file = file.path(output_dir, "pca_coordinates_3patients.csv"),
  row.names = FALSE
)

# PCA plot
pca_plot <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = condition, shape = patient)
) +
  geom_point(size = 4) +
  geom_text_repel(
    aes(label = sample),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.4,
    min.segment.length = 0
  ) +
  xlab(paste0("PC1 (", pc1_var, "% variance)")) +
  ylab(paste0("PC2 (", pc2_var, "% variance)")) +
  ggtitle("PCA of RNA-seq samples") +
  labs(
    color = "Condition",
    shape = "Patient ID"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(15, 30, 15, 30)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(output_dir, "pca_3patients.png"),
  plot = pca_plot,
  width = 7,
  height= 5,
  dpi = 300
)

# calculates sample to sample distances
sample_dists <- dist(t(vsd_mat))
sample_dist_matrix <- as.matrix(sample_dists)

write.csv(
  x = sample_dist_matrix,
  file = file.path(output_dir, "sample_distance_matrix_3patients.csv")
)

# sample dist. heatmap
png(
  filename = file.path(output_dir, "sample_distance_heatmap_3patients.png"),
  width = 1800,
  height = 1600,
  res = 300
)

pheatmap(
  mat = sample_dist_matrix,
  main = "Sample to sample distances",
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists
)

dev.off()

# MA plot 
# MA plot data
res <- results(
  dds,
  contrast = c("condition", "tumor", "normal")
)

ma_df <- as.data.frame(res)
ma_df$gene_id <- rownames(ma_df)

ma_df <- ma_df[
  !is.na(ma_df$baseMean) &
    !is.na(ma_df$log2FoldChange) &
    !is.na(ma_df$padj),
]

# MA plot groups
ma_df$ma_group <- "not_significant"
ma_df$ma_group[ma_df$padj < 0.05 & ma_df$log2FoldChange > 0] <- "up_in_tumor"
ma_df$ma_group[ma_df$padj < 0.05 & ma_df$log2FoldChange < 0] <- "down_in_tumor"

ma_plot <- ggplot(
  ma_df,
  aes(x = baseMean, y = log2FoldChange, color = ma_group)
) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_x_log10() +
  coord_cartesian(ylim = c(-6, 6)) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "not_significant" = "grey70",
      "up_in_tumor" = "blue",
      "down_in_tumor" = "red"
    ),
    labels = c(
      "not_significant" = "Not significant",
      "up_in_tumor" = "Upregulated in tumor",
      "down_in_tumor" = "Downregulated in tumor"
    )
  ) +
  labs(
    title = "MA plot: tumor vs normal",
    x = "Mean of normalized counts",
    y = "Log2 fold change",
    color = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(output_dir, "ma_plot_3patients_custom.png"),
  plot = ma_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# summary
cat("script completed\n")
cat("Number of samples:", ncol(vsd_mat), "\n")
cat("Number of genes:", nrow(vsd_mat), "\n")
cat("Results written to:", output_dir, "\n")
