library(rtracklayer)
gtf_file <- "reference/gencode_v49_salmon/gencode.v49.annotation.gtf.gz"
output_file <- "reference/gencode_v49_salmon/tx2gene.tsv"

gtf <- import(gtf_file)

gtf_table <- as.data.frame(gtf)

transcript_rows <- gtf_table[gtf_table$type == "transcript",]

tx2gene <-transcript_rows[, c("transcript_id", "gene_id")]

tx2gene <- unique(tx2gene)

write.table(
  tx2gene,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)