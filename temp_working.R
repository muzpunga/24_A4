library(DESeq2)
library(ggplot2)
library(apeglm)

# ============================================================================
# Create figures directory if it does not already exist
# ============================================================================

dir.create("figures", showWarnings = FALSE)

# ============================================================================
# Read six count files into one count matrix
# ============================================================================

condition_a_files <- c(
  "sample_024_counts.txt",
  "sample_026_counts.txt",
  "sample_029_counts.txt"
)

condition_b_files <- c(
  "sample_007_counts.txt",
  "sample_017_counts.txt",
  "sample_032_counts.txt"
)

counts_dir <- "../../shared/SLE777/A4_counts"

dir(counts_dir)

condition_a_files %in% dir(counts_dir)
condition_b_files %in% dir(counts_dir)

all_files <- c(condition_a_files, condition_b_files)

available_files <- intersect(
  all_files,
  dir(counts_dir)
)

available_files

# ============================================================================
# Build sample metadata
# ============================================================================

sample_metadata <- data.frame(
  sample_id = available_files,
  condition = ifelse(
    available_files %in% condition_a_files,
    "A",
    "B"
  ),
  stringsAsFactors = FALSE
)

# Set condition A as reference level

sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("A", "B")
)

rownames(sample_metadata) <- sample_metadata$sample_id

sample_metadata

# ============================================================================
# Read count files
# ============================================================================

count_list <- lapply(
  available_files,
  function(f) {
    
    read.table(
      file.path(counts_dir, f),
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  }
)

# ============================================================================
# Build count matrix
# ============================================================================

count_matrix <- data.frame(
  gene_id = count_list[[1]][, 1]
)

for (i in seq_along(count_list)) {
  
  count_matrix[[available_files[i]]] <- count_list[[i]][, 2]
  
}

rownames(count_matrix) <- count_matrix$gene_id
count_matrix$gene_id <- NULL

count_matrix <- as.matrix(count_matrix)

dim(count_matrix)

# ============================================================================
# Filter low-count genes
#
# Keep genes with CPM > 1 in at least 3 samples.
# Three samples is the size of the smallest experimental group.
# ============================================================================

library_sizes <- colSums(count_matrix)

round(
  library_sizes / 1e6,
  2
)

cpm <- count_matrix

for (i in seq_len(ncol(count_matrix))) {
  
  cpm[, i] <-
    count_matrix[, i] /
    library_sizes[i] *
    1e6
  
}

# CPM should sum to approximately one million per sample
# colSums(cpm)

keep <- rowSums(cpm > 1) >= 3

table(keep)

counts_filtered <- count_matrix[keep, ]
cpm_filtered <- cpm[keep, ]

dim(count_matrix)
dim(counts_filtered)

# ============================================================================
# Compute log-CPM
#
# +1 is a standard pseudocount used to avoid log(0)
# ============================================================================

log_cpm <- log2(cpm_filtered + 1)

dim(log_cpm)

head(log_cpm)

# ============================================================================
# Boxplots before and after normalisation
# ============================================================================

boxplot(
  log2(counts_filtered + 1),
  las = 2,
  ylab = "log2(count + 1)",
  main = "Before normalisation"
)

boxplot(
  log_cpm,
  las = 2,
  ylab = "log2(CPM + 1)",
  main = "After normalisation"
)

# ============================================================================
# Principal Component Analysis (PCA)
# ============================================================================

dim(log_cpm)

# PCA requires samples in rows and genes in columns

dim(t(log_cpm))

pca <- prcomp(t(log_cpm))

dim(pca$x)

head(rownames(pca$x))

# ============================================================================
# Percentage variance explained
#
# percent_var[1] is the percentage variance explained by PC1
# ============================================================================

percent_var <- round(
  100 * pca$sdev^2 /
    sum(pca$sdev^2),
  1
)

percent_var[1:4]

# ============================================================================
# PCA plotting data frame
# ============================================================================

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

pca_df <- merge(
  pca_df,
  sample_metadata,
  by = "sample_id"
)

# ============================================================================
# PCA plot
# ============================================================================

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = condition
  )
) +
  geom_point(size = 4) +
  geom_text(
    aes(label = sample_id),
    vjust = -1.2,
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "A" = "#66c2a5",
      "B" = "#fc8d62"
    )
  ) +
  labs(
    title = "PCA of normalised expression",
    x = paste0(
      "PC1 (",
      percent_var[1],
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      percent_var[2],
      "% variance)"
    ),
    colour = "Condition"
  ) +
  theme_bw()

p_pca

# Save PCA

ggsave(
  filename = "figures/pca_plot.png",
  plot = p_pca,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# ============================================================================
# DESeq2 analysis
# ============================================================================

all(
  colnames(counts_filtered) ==
    sample_metadata$sample_id
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_metadata,
  design = ~ condition
)

dds <- DESeq(dds)

class(dds)

levels(dds$condition)

resultsNames(dds)

# ============================================================================
# Extract results for condition B versus condition A
#
# Positive log2FC = higher in condition B
# Negative log2FC = higher in condition A
# ============================================================================

res <- results(
  dds,
  contrast = c(
    "condition",
    "B",
    "A"
  )
)

summary(res)

# ============================================================================
# Library sizes and size factors
# ============================================================================

sizeFactors(dds)

library_sizes <- colSums(counts_filtered)

round(
  library_sizes /
    mean(library_sizes),
  3
)

round(
  sizeFactors(dds),
  3
)

# ============================================================================
# Dispersion plot
# ============================================================================

plotDispEsts(dds)

png(
  filename = "figures/dispersion_plot.png",
  width = 7,
  height = 6,
  units = "in",
  res = 600
)

plotDispEsts(dds)

dev.off()

# ============================================================================
# Shrink fold changes
#
# Shrinkage reduces unstable fold changes for low-count genes
# while preserving strong, well-supported effects.
# ============================================================================

res_shrunk <- lfcShrink(
  dds,
  coef = "condition_B_vs_A",
  type = "apeglm"
)

summary(res_shrunk)

# Convert to data frame for plotting

res_df <- as.data.frame(res_shrunk)

# ============================================================================
# Volcano plot preparation
# ============================================================================

res_df$neg_log10_padj <- -log10(res_df$padj)

# ============================================================================
# Classify significant genes
#
# Significant:
#   padj < 0.05
#   |log2FC| > 1
#
# Positive log2FC = Up in B
# Negative log2FC = Up in A
# ============================================================================

res_df$status <- "Not significant"

res_df$status[
  res_df$padj < 0.05 &
    res_df$log2FoldChange > 1
] <- "Up in B"

res_df$status[
  res_df$padj < 0.05 &
    res_df$log2FoldChange < -1
] <- "Up in A"

table(res_df$status)

# ============================================================================
# Volcano plot
# ============================================================================

p_volcano <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    colour = status
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  scale_colour_manual(
    values = c(
      "Up in B" = "#d73027",
      "Up in A" = "#4575b4",
      "Not significant" = "grey70"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    colour = "grey40"
  ) +
  labs(
    title = "Differential expression: B versus A",
    x = "log2 fold change (B / A)",
    y = "-log10 adjusted p-value",
    colour = ""
  ) +
  theme_bw()

p_volcano

# Save volcano plot

ggsave(
  filename = "figures/volcano_plot.png",
  plot = p_volcano,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)
