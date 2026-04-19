# libs
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

# file paths
all_results_file <- "results/deseq2_results_summary_3patients/all_results_cleaned_3patients.csv"
up_file <- "results/deseq2_results_summary_3patients/upregulated_in_tumor_3patients.csv"
down_file <- "results/deseq2_results_summary_3patients/downregulated_in_tumor_3patients.csv"
output_dir <- "results/enrichment_3patients"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# loads tables
all_results <- read.csv(all_results_file, stringsAsFactors = FALSE)
up_results <- read.csv(up_file, stringsAsFactors = FALSE)
down_results <- read.csv(down_file, stringsAsFactors = FALSE)

# keeps only tested genes
tested_results <- all_results[!is.na(all_results$padj),, drop = FALSE]

# checks for required columns (issues with ensembl id's and gene symbols)
required_columns <- c("ensembl_gene_id_clean", "gene_symbol")
for (col_name in required_columns) {
  if (!(col_name %in% colnames(tested_results))) {
    stop(paste("missing column in tested results: ", col_name))
  } 
  if (!(col_name %in% colnames(up_results))) {
    stop(paste("missing column in upregulated results: ", col_name))
  }
  if (!(col_name %in% colnames(down_results))) {
    stop(paste("missing column in downregulated results:", col_name))
  }
}

# maps cleaned Ensembl ids to Entrez ids
map_to_entrez <- function(df, group_name) {
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = df$ensembl_gene_id_clean,
    column  = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  mapped_df <- data.frame(
    ensembl_gene_id_clean = df$ensembl_gene_id_clean,
    gene_symbol = df$gene_symbol,
    entrez_id = unname(entrez_ids),
    group = group_name,
    stringsAsFactors = FALSE
  )
  
  mapped_df <- mapped_df[!is.na(mapped_df$entrez_id),, drop = FALSE]
  mapped_df <- mapped_df[mapped_df$entrez_id != "",, drop = FALSE]
  mapped_df <- mapped_df[!duplicated(mapped_df$entrez_id), ,drop = FALSE]
  
  mapped_df
}

# makes mapped gene tables
universe_map <- map_to_entrez(tested_results, "tested_background")
up_map <- map_to_entrez(up_results, "up_in_tumor")
down_map <- map_to_entrez(down_results, "down_in_tumor")

universe_entrez <- universe_map$entrez_id
up_entrez <- up_map$entrez_id
down_entrez <- down_map$entrez_id

# saves mapped gene list
write.csv(
  universe_map,
  file = file.path(output_dir, "tested_background_mapped_genes.csv"),
  row.names = FALSE
)

write.csv(
  up_map,
  file = file.path(output_dir, "upregulated_mapped_genes.csv"),
  row.names = FALSE
)

write.csv(
  down_map,
  file = file.path(output_dir, "downregulated_mapped_genes.csv"),
  row.names = FALSE
)

# runs ORA
run_enrichment <- function(gene_ids, universe_ids) {
  go_bp <- enrichGO(
    gene = gene_ids,
    universe = universe_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  reactome <- enrichPathway(
    gene = gene_ids,
    universe = universe_ids,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  list(
    go_bp = go_bp,
    reactome = reactome
  )
}

up_enrichment <- run_enrichment(up_entrez, universe_entrez)
down_enrichment <- run_enrichment(down_entrez, universe_entrez)

# saves enrichment tables
save_enrichment_table <- function(enrichment_object, out_file) {
  enrichment_table <- as.data.frame(enrichment_object)
  
  if (nrow(enrichment_table) == 0) {
    enrichment_table <- data.frame(
      note = "no siginificant enrichment terms found",
      stringsAsFactors = FALSE
    )
  }
  
  write.csv(
    enrichment_table,
    file = out_file,
    row.names = FALSE
  )
}

save_enrichment_table(
  up_enrichment$go_bp,
  file.path(output_dir, "go_bp_upregulated_in_tumor.csv")
)

save_enrichment_table(
  down_enrichment$go_bp,
  file.path(output_dir, "go_bp_downregulated_in_tumor.csv")
)

save_enrichment_table(
  up_enrichment$reactome,
  file.path(output_dir, "reactome_upregulated_in_tumor.csv")
)

save_enrichment_table(
  down_enrichment$reactome,
  file.path(output_dir, "reactome_downregulated_in_tumor.csv")
)

# saves enrichment objects
saveRDS(
  up_enrichment$go_bp,
  file = file.path(output_dir, "go_bp_upregulated_in_tumor.rds")
)

saveRDS(
  down_enrichment$go_bp,
  file = file.path(output_dir, "go_bp_downregulated_in_tumor.rds")
)

saveRDS(
  up_enrichment$reactome,
  file = file.path(output_dir, "reactome_upregulated_in_tumor.rds")
)

saveRDS(
  down_enrichment$reactome,
  file = file.path(output_dir, "reactome_downregulated_in_tumor.rds")
)

# helper to save dotplots
save_dotplot <- function(enrichment_object, plot_title, out_file) {
  enrichment_table <- as.data.frame(enrichment_object)
  
  if (nrow(enrichment_table) == 0) {
    return(NULL)
  }
  
  plot_object <- dotplot(enrichment_object, showCategory = 15) +
    ggtitle(plot_title)
  
  ggsave(
    filename = out_file,
    plot = plot_object,
    width = 10,
    height = 14,
    dpi = 300
  )
}

# plots
save_dotplot(
  up_enrichment$go_bp,
  "GO Biological Process enrichment: upregulated in tumor",
  file.path(output_dir, "go_bp_upregulated_in_tumor_dotplot.png")
)

save_dotplot(
  down_enrichment$go_bp,
  "GO Biological Process enrichment: downregulated in tumor",
  file.path(output_dir, "go_bp_downregulated_in_tumor_dotplot.png")
)

save_dotplot(
  up_enrichment$reactome,
  "Reactome enrichment: upregulated in tumor",
  file.path(output_dir, "reactome_upregulated_in_tumor_dotplot.png")
)

save_dotplot(
  down_enrichment$reactome,
  "Reactome enrichment: downregulated in tumor",
  file.path(output_dir, "reactome_downregulated_in_tumor_dotplot.png")
)

#summary
count_terms <- function(enrichment_object) {
  nrow(as.data.frame(enrichment_object))
}

summary_lines <- c(
  "Enrichment analysis summary",
  paste("Tested background genes mapped to Entrez:", length(universe_entrez)),
  paste("Upregulated genes mapped to Entrez: " ,length(up_entrez)),
  paste("Downregulated genes mapped to Entrez:", length(down_entrez)),
  paste("GO BP terms enriched in upregulated genes:", count_terms(up_enrichment$go_bp)),
  paste("GO BP terms enriched in downregulated genes:", count_terms(down_enrichment$go_bp)),
  paste("Reactome terms enriched in upregulated genes:", count_terms(up_enrichment$reactome)),
  paste("Reactome terms enriched in downregulated genes:", count_terms(down_enrichment$reactome))
)

writeLines(
  summary_lines,
  con = file.path(output_dir, "enrichment_summary_3patients.txt")
)

summary_lines
