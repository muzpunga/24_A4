# All code chunks in 24_report.Rmd must be run prior.

group_id         <- params$group_id
gene_of_interest <- params$gene_of_interest
n_sig            <- nrow(sig_genes)
n_up_in_B        <- nrow(up_in_b) # R is case sensitive
n_down_in_B      <- nrow(down_in_b)
top_gene_id      <- top_gene
goi_log2fc       <- goi_lfc
goi_padj         <- goi_padj
goi_direction    <- strsplit(tolower(goi_status), " ")[[1]][1]
ego_df           <- as.data.frame(ego) # needed for next line
top_go_bp_id     <- ego_df$ID[which.min(ego_df$p.adjust)]

results_csv <- read.csv("24_results.csv", stringsAsFactors = FALSE)

results_csv[, 2] <- c(
  group_id,
  gene_of_interest,
  n_sig,
  n_up_in_B,
  n_down_in_B,
  top_gene_id,
  goi_log2fc,
  goi_padj,
  goi_direction,
  top_go_bp_id
)

write.csv(results_csv, "24_results.csv", row.names = FALSE)