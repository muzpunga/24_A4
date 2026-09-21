library(DESeq2)
library(ggplot2)

# Read six count files into one count matrix
condition_a_files <- c("sample_024_counts.txt", "sample_026_counts.txt", "sample_029_counts.txt")
condition_b_files <- c("sample_007_counts.txt", "sample_017_counts.txt", "sample_032_counts.txt")

counts_dir <- "../../shared/SLE777/A4_counts"

dir(counts_dir)

condition_a_files %in% dir(counts_dir)
condition_b_files %in% dir(counts_dir)


all_files <- c(condition_a_files, condition_b_files)
available_files <- intersect(all_files, dir(counts_dir))
available_files

# Build sample metadata
sample_metadata <- data.frame(
  sample_id = available_files,
  condition = ifelse(available_files %in% condition_a_files, "A", "B"),
  stringsAsFactors = FALSE)

# Build sample metadata with A as reference level
sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("A", "B"))

rownames(sample_metadata) <- sample_metadata$sample_id

count_list <- lapply(available_files, function(f) {
    dat <- read.table(
      file.path(counts_dir, f),
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE)
    
    dat
  }
)

# Read six count files into one count matrix
count_matrix <- data.frame(
  gene_id = count_list[[1]][, 1]
)

for(i in seq_along(count_list)) {
  count_matrix[[available_files[i]]] <- count_list[[i]][, 2]
}

rownames(count_matrix) <- count_matrix$gene_id
count_matrix$gene_id <- NULL

count_matrix <- as.matrix(count_matrix)

# A widely used rule is to keep a gene if it has a CPM above 1 in at least as many samples as the smallest group
# Week 8 6. Filtering low-count genes


library_sizes <- colSums(count_matrix)
round(library_sizes / 1e6, 2)

cpm <- count_matrix

for (i in 1:ncol(count_matrix)) {
  cpm[, i] <- count_matrix[, i] / library_sizes[i] * 1e6
}

# must sum to 1 million
# colSums(cpm)

# Filter unexpressed genes 1.
keep <- rowSums(cpm > 1) >= 3
table(keep)
# Filter unexpressed genes 2.
counts_filtered <- count_matrix[keep, ] 
cpm_filtered <- cpm[keep, ]
dim(count_matrix)
dim(counts_filtered)

# Compute log-CPM 1. +1 is standard approach to handling zeros
log_cpm <- log2(cpm_filtered + 1)
dim(log_cpm)
head(log_cpm)

boxplot(
  log2(counts_filtered + 1),
  las = 2,
  ylab = "log2(count + 1)",
  main = "Before normalisation")

boxplot(
  log_cpm,
  las = 2,
  ylab = "log2(CPM + 1)",
  main = "After normalisation")



# PCA 1. checks
dim(log_cpm)
# transposition is required - after transposition the row count should be small and colum count high
dim(t(log_cpm))
# PCA 2. perform PCA
pca <- prcomp(t(log_cpm))
dim(pca$x)
head(rownames(pca$x))


# Report % variance explained by PC1
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
percent_var[1:4] # percent_var[1] is the variance explained by PC1

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2])

pca_df <- merge(pca_df, sample_metadata, by = "sample_id")


p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 4) +
  scale_colour_manual(values = c("A" = "#66c2a5", "B" = "#fc8d62")) +
  labs(
    title = "PCA of normalised expression",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    colour = "Stage"
  ) + theme_bw() +
  geom_text(aes(label = sample_id), vjust = -1.2, size = 3, show.legend = FALSE)

# Build DESeqDataSet with ~ condition 1. checks
all(colnames(counts_filtered) == sample_metadata$sample_id) #stopifnot()
# Build DESeqDataSet with ~ condition 2. create SESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_metadata,
  design = ~ condition)
# Build DESeqDataSet with ~ condition 3. convert to dds and check its object type
dds <- DESeq(dds)
class(dds)

levels(dds$condition)
resultsNames(dds) # needed for Shrink fold below


# extract the results for the contrast of condition B versus condition A
res <- results(
  dds,
  contrast = c("condition", "B", "A"))

summary(res)
# Differential expression analysis was performed using DESeq2 with condition A specified as the reference level (levels = c("A", "B")). Results were extracted for the contrast B versus A (contrast = c("condition", "B", "A")). Consequently, the reported log2 fold changes represent expression in condition B relative to condition A. A positive log2 fold change indicates higher expression in condition B, whereas a negative log2 fold change indicates higher expression in condition A.


dds_esf <- estimateSizeFactors(dds)
sizeFactors(dds_esf)
# Condition A appears to be under represented, and B over-represented

library_sizes <- colSums(counts_filtered)
round(library_sizes / mean(library_sizes), 3)
round(sizeFactors(dds_esf), 3)

# Dispersion

dds_disp <- estimateDispersions(dds_esf)
plotDispEsts(dds_disp)
# The dispersion plot showed the expected inverse relationship between mean normalized counts and dispersion. Gene-wise dispersion estimates (black points) followed the fitted dispersion trend (red curve), while empirical Bayes shrinkage produced final dispersion estimates (blue points) that were stabilized towards the fitted trend. Overall, the dispersion model appeared well-behaved with no evidence of major fitting problems.


# Shrink fold changes using coefficient from resultsNames()
library(apeglm)
res_shrunk <- lfcShrink(
  dds,
  coef = "condition_B_vs_A", # the coefficient name that resultsNames() reports
  type = "apeglm")
summary(res_shrunk)

head(as.data.frame(res_shrunk))

