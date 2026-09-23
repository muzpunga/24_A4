####################################
# Working analysis script for assessment 4
# Samples: a (10 days) vs b (45 days), flor yeast
# Author: Ivy Craw
# Date: 2026-09-21
####################################

# In our Bash pipeline, we used the command line to process RNA-seq data
# We downloaded the raw FASTQ files, checked their quality, aligned reads to a reference genome
# ...converted alignment files into BAM format, sorted and indexed them
# ...and finally counted how many reads mapped to each gene
# The output of the read-counting step is a gene count file, containing gene ID and n of reads assigned to that gene in one sample per row

# For this assessment, we are using pre-generated count files assigned to our group

# ======= Step 0 Libraries =======

# Load ggplot2 package for plotting
library(ggplot2)
# Load ggpel for dynamic labelling, okayed by unit chair for use OPTIONAL
library(ggrepel)
# Load DESeq2, DESeq2 comes from Bioconductor rather than CRAN
library(DESeq2)
# For final ...
library(pheatmap)
# These are used for...
library(org.Sc.sgd.db)
# Used for...
library(clusterProfiler)

# ======= Step 0.5 Set directory structure =======
# Already set as Git repo
# Enter the main project directory
cd 24_A4_working

####################################
# Build the count matrix
####################################

# Combine six HTSeq-count files into one gene x sample matrix.
# Samples: condition a versus condition b, for assessment purposes

# To compare gene expression between conditions, we need all six samples side by side in one table
# The end product of this week is a count matrix
# One row per gene and one column per sample and a read count in every cell
# We also build a sample metadata table
# These two objects are required for differential expression analysis (DESeq2)

### Descriptions of samples
# Unknown. Counts generated for assessment purposes. 

####################################
# The sample metadata table
####################################

# Create sample metadata table with sample identifiers and condition
sample_metadata <- data.frame(
  sample_id = c("sample_024", "sample_026", "sample_029",
                "sample_007", "sample_017", "sample_032"),
  condition = c("A", "A", "A",
                "B", "B", "B")
)
sample_metadata

### Conditions should be factors
# Statistical models in R need categorical variables stored as factors
# At the moment condition is a character column
class(sample_metadata$condition)
str(sample_metadata)

# We want condition a to be the reference
# ...so that later a positive fold change will mean "higher in b than a"

# Convert the condition column to a factor with a as the reference level
sample_metadata$condition <- factor(
  sample_metadata$condition,
  levels = c("A", "B")
)
# Check that the conversion worked
levels(sample_metadata$condition)
# Check the structure of the metadata table
str(sample_metadata)

### Row names

# Name rows of the metadata after the samples
rownames(sample_metadata) <- sample_metadata$sample_id

# Confirm the metadata table is correct after these changes
sample_metadata

### Finding the count files

# Locate the count files assigned for the assessment
# We need to be specific here as there are many files
condition_A_files <- c(
  "sample_024_counts.txt",
  "sample_026_counts.txt",
  "sample_029_counts.txt"
)

condition_B_files <- c(
  "sample_007_counts.txt",
  "sample_017_counts.txt",
  "sample_032_counts.txt"
)

# Define the file path for our files
# This is the only time we are using an absolute path
# This was called for in the assessment instructions
counts_dir <- "/home/shared/SLE777/A4_counts"

# Create vectors to save the files and their paths
count_files <- c(condition_A_files, condition_B_files)
count_paths <- file.path(counts_dir, count_files)

# Confirm that you can see the full file paths and file names
count_paths

# Safety checks to stop if not all six samples
stopifnot(length(count_files) == 6)
stopifnot(all(file.exists(count_paths)))

### Deriving sample names from file names

# Check file names again for what to strip
count_files

# Strip until we have sample identifiers only
sample_ids <- gsub("_counts\\.txt$", "", count_files)
sample_ids

# Check that every sample id derived from file name also appears in metadata table
# Stop the script if all results do not evaluate as TRUE
stopifnot(all(sample_ids %in% sample_metadata$sample_id))

####################################
# Merging the count files into one table
####################################

# Each count file is a two-column table: 
# ...a gene identifier, and the count for that one sample 
# To put six samples side by side
# ...we need to join these tables so that each row still refers to the same gene 

### Reading all six files with a loop
# The if else was added to the loop to create the counts_table object from scratch
# This deviates from the practical code as that included demonstration code
for (i in 1:length(count_paths)) {
  
  one_sample <- read.table(
    count_paths[i],
    header = FALSE,
    col.names = c("gene_id", sample_ids[i])
  )
  
  # Quick checks to reflect structure set earlier in code
  # Stop if the count table being processed does not have two columns
  # or if the first column is not character
  # or if the second column is not numeric
  stopifnot(
    ncol(one_sample) == 2,
    is.character(one_sample$gene_id),
    is.numeric(one_sample[[2]])
  )

  # Merge the new sample into the growing counts table
  # The first time through the loop the counts_table does not yet exist
  if (i == 1) {
    counts_table <- one_sample
  } else {
    counts_table <- merge(counts_table, one_sample, by = "gene_id")
  }
}

### Quick checks and guardrails for counts_table creation

# Confirm merged table contains one gene_id column plus six sample columns
stopifnot(ncol(counts_table) == length(sample_ids) + 1)

# Column names should be gene ID plus our sample IDs
expected_colnames <- c("gene_id", sample_ids)
stopifnot(all(colnames(counts_table) == expected_colnames))

# Check that the final table to see how many gene IDs and samples you have
# Remember that these are still inclusive of summary rows
message(
  "Merged count table created successfully: ",
  nrow(counts_table), " rows x ",
  ncol(counts_table), " columns."
)

####################################
# Separating the HTSeq-count summary rows
####################################

# Every count file has summary rows
# These were merged along with gene counts and will need to be separated

# Apply grepl to the wider table
is_summary <- grepl("^__", counts_table$gene_id)
sum(is_summary)

# Split table in two
htseq_summary <- counts_table[is_summary, ]
gene_table <- counts_table[!is_summary, ]

# Check that the split worked
dim(counts_table)

# Check that summary rows are in the summary table
dim(htseq_summary)

# Check that the gene table has no summary rows
dim(gene_table)

### Guardrails to be confident we have removed them

# Check that the summary table has five rows, one for each of the five HTSeq-count summary categories
stopifnot(nrow(htseq_summary) == 5)

# Check that the gene table has no summary rows
stopifnot(!any(grepl("^__", gene_table$gene_id)))

# And lastly, check that the gene table has the same number of rows as the original count files
# ...because we have not lost any genes in the merge
# The number of rows in the gene table should match the number of rows in any one of the original count files, minus the five summary rows
stopifnot(nrow(gene_table) == nrow(one_sample) - 5)

### These rows tell you how many reads were discarded because they fell outside any gene
# ...because they were ambiguous between overlapping genes
# or because their alignment quality was too low
# A sample with an unusually high __no_feature count deserves a closer look
htseq_summary

####################################
# Building the count matrix
####################################

# Need to remove gene_id so there is no text in table
# This is so that we can convert the table to a numeric matrix
count_matrix <- as.matrix(gene_table[, -1])
rownames(count_matrix) <- gene_table$gene_id

# Manually check the new object has the same dimensions minus one column
# Inspect visually the structure and first few rows
dim(count_matrix)
str(count_matrix)
head(count_matrix)

# And guardrailds for safety
stopifnot(
  nrow(count_matrix) == nrow(gene_table),
  ncol(count_matrix) == ncol(gene_table) - 1,
  all(rownames(count_matrix) == gene_table$gene_id),
  all(colnames(count_matrix) == sample_ids),
  is.matrix(count_matrix),
  is.numeric(count_matrix)
)

####################################
# Checking that the samples line up
####################################

# The columns of the count matrix and the rows of the metadata 
# ...must describe the same samples in the same order
# Nothing enforces this so we watch for silent errors

# First look at the two sets of names side by side
colnames(count_matrix)
sample_metadata$sample_id

# Now force the columns into the metadata order and check
count_matrix <- count_matrix[, sample_metadata$sample_id]
all(colnames(count_matrix) == sample_metadata$sample_id)
identical(colnames(count_matrix), rownames(sample_metadata))


# Also get R to check for us
stopifnot(all(colnames(count_matrix) == sample_metadata$sample_id))
stopifnot(identical(colnames(count_matrix), rownames(sample_metadata)))

### Always reorder the matrix columns to match the metadata rows
# ...always assert that they match
# ...and never continue past a failed check


####################################
# Library sizes and read assignments
####################################

# The library size of a sample is the total number of reads assigned to genes
# Samples with more reads will show higher counts for every gene
# ...which is why we cannot compare raw counts between samples
# Normalising converts counts into a quantity that does not depend on how deeply the sample was sequenced

# Remember to use colSums
# A lib size is a total per sample, and samples are columns in our matrix
library_sizes <- colSums(count_matrix)
length(library_sizes)
library_sizes

# Scale library sizes to millions of reads and round to 2 decimal places
round(library_sizes / 1e6, 2)

### Summary of reads that were not assigned to any gene

# Summary rows we set aside earlier
# Convert to matrix, removing header, then add up each column
summary_matrix <- as.matrix(htseq_summary[, -1])
rownames(summary_matrix) <- htseq_summary$gene_id
summary_matrix

# Sum the rows of the summary matrix to get the total number of reads that were not assigned to any gene
unassigned <- colSums(summary_matrix)
unassigned

### Calculate the proportion of reads assigned to genes

# First build a data frame
qc_table <- data.frame(
  sample_id = names(library_sizes),
  library_size = as.numeric(library_sizes),
  unassigned = as.numeric(unassigned[names(library_sizes)])
)
qc_table

# Add two more columns using existing columns
qc_table$total_reads <- qc_table$library_size + qc_table$unassigned
qc_table$percent_assigned <- 100 * qc_table$library_size / qc_table$total_reads
qc_table

# A healthy RNA-seq sample from a well-annotated genome such as yeast should assign a large majority of its reads to genes
# A sample that is markedly worse than the others is a candidate for exclusion
# Should a least be flagged in report

### Plotting the library sizes

# We need to merge condition alongside the library size for coloured bars
qc_table <- merge(qc_table, sample_metadata, by = "sample_id")
qc_table

# Fill aesthetics maps a column to bar colour
# theme_bw() function replaces the grey default background
# theme() call rotates the x-axis labels so that the long sample names do not overlap
ggplot(qc_table, aes(x = sample_id, y = library_size / 1e6, fill = condition)) +
  geom_col() +
  labs(
    title = "Library size per sample",
    x = "Sample",
    y = "Assigned reads (millions)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Plot percent_assigned for each sample as a bar chart, coloured by condition, with a sensible title and axis labels
ggplot(qc_table, aes(x = sample_id, y = percent_assigned, fill = condition)) +
  geom_col() +
  labs(
    title = "Percentage of reads assigned to genes",
    x = "Sample",
    y = "Reads assigned (%)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####################################
# A first look across samples
####################################

# RNA-seq counts often have a skewed distribution
# Many genes have low counts while a smaller number have very high counts
# Log transformation makes the distributions easier to compare across samples

# Draw histogram of counts across samples
# The las = 2 argument turns the sample labels sideways so they fit
# Logarithm cannot take zero so we add 1 across whole matrix
boxplot(
  log10(count_matrix + 1),
  las = 2,
  ylab = "log10(count + 1)",
  main = "Distribution of gene counts per sample"
)

####################################
# Filtering genes with no expression
####################################

# Many genes are not expressed at all in this experiment
# They carry no information, they slow the analysis down
# ...and they cost us statistical power when we correct for multiple testing in Week 9

### Count how many genes have zero reads across every sample

# Create new var for total n of rows count as total number of genes
gene_totals <- rowSums(count_matrix)

### Remove the genes with zero reads across samples

# Create new variable to store filtered new matrix
# Note, a logical vector goes before the comma in the filter
# blanks keeps all the columns
count_matrix_filtered <- count_matrix[gene_totals > 0, ]

# Calculate summary statistics of our genes across samples
dim(count_matrix)
# Summary statistics for only genes with reads
dim(count_matrix_filtered)

# Quick activity to record genes with at least ten reads in total across samples
sum(gene_totals >= 10)
count_matrix_10 <- count_matrix[gene_totals >= 10, ]
dim(count_matrix_10)

####################################
# Week 7 summary
####################################

# We checked what our samples actually were
# Joined six count files into one table with merge()
# Separated the summary rows
# Built a numeric count matrix
# Checked that the samples lined up 
# ...and looked at library sizes and assignment rates

####################################
# END OF FIRST SECTION             #
####################################


##############################################
# Normalisation and exploring the experiment
##############################################

# Define our library sizes again
library_sizes <- colSums(count_matrix)
round(library_sizes / 1e6, 2)

### 5. Counts per million

# divide each count by its sample's library size
# then multiply by one million so the numbers are readable 
# We now do it to every cell

# Each sample has its own library size
# Work through columns one at a time

# Copy count matrix so dimensions are the same
cpm <- count_matrix

# For each column i take col of counts
# divide it by that sample's lib size
# multiply by a million
# and put answer back in the count matrix copy cpm
for (i in 1:ncol(count_matrix)) {
  cpm[, i] <- count_matrix[, i] / library_sizes[i] * 1e6
}

# Check that we still have x genes and x samples
dim(cpm)
# Inspect results
head(cpm)
# and confirm that every column of CPM matrix sum is exactly one mill 1e+06
colSums(cpm)

# CPM corrects for sequencing depth
# It does not correct for gene length 
# ...so CPM values cannot be compared between different genes
# They can be compared for the same gene across samples
# ...which is exactly what differential expression needs

### 6. Filtering low-count genes

# We cannot reliably test genes with almost no reads anywhere
# They inflate the multiple-testing correction we apply later
# costing us the power to detect genes of interest

# A widely used rule is to keep a gene if it has a CPM above 1 in at least as many samples as the smallest group
# Our smallest group has three replicates

keep <- rowSums(cpm > 1) >= 3
table(keep)

# The expression cpm > 1 gives a boolean matrix the same shape as cpm
# R treats TRUE as 1, so rowSums() of that matrix counts, for each gene, how many samples passed
# We then keep the genes where that count is at least three

# Apply the filter to both matrices
counts_filtered <- count_matrix[keep, ]
cpm_filtered <- cpm[keep, ]
dim(count_matrix)
dim(counts_filtered)

# Filtering through a logical vector before the comma, all columns kept
# Note that we keep the filtered raw counts as well
# DESeq2 next week wants raw counts, not CPM, because it models the counting process itself
# The CPM matrix is for exploration and plotting only

### 7. Log transformation

# Count data are strongly skewed
# A handful of genes carry enormous counts and dominate every plot
# Taking logarithms compresses that range
# As always, we add 1 first, because the log of zero is undefined

log_cpm <- log2(cpm_filtered + 1)

# Inspect the log CPM matrix
dim(log_cpm)
head(log_cpm)

# Earlier in the srcipt we used log10 to squash the y-axis of a boxplot
# Here we switch to log2, and the change is deliberate
# Base 2 is conventional in genomics because a difference of 1 on the log2 scale means a doubling of expression
# ...a difference of 2 means a fourfold change, and so on
# That interpretation carries directly into the log2 fold changes we calculate later in the analysis
# For simply squashing a plot's axis either base works; for talking about fold changes, base 2 is the one that lets you read the number

### 8. Did the normalisation work?

# We have claimed that normalisation makes samples comparable
# Do not take it on trust. Draw the distributions before and after and look.
# Because both matrices are numeric, base R gives one box per sample from a single line of code

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

##################################
# Principal component analysis ###
##################################

# Let's consider
# if the two stages really differ
# ...samples from the same stage should look more like each other than like samples from the other stage

# Our data have thousands of dimensions, one per gene, and we cannot draw that
# PCA finds a small number of new axes, called principal components, that capture as much of the variation as possible
# It lets us show six samples measured on thousands of genes on an ordinary scatter plot
# PCA needs parts per million

# Each component is built from a combination of many genes
# PC1 is the direction along which the samples differ most
# PC2 the next most, and so on
# If the experiment worked, PC1 usually separates the conditions

# 9.1 Transposing
# Risk of silent error here! 

# Transpose log-CPM matrix so samples are rows and genes are columns for prcomp()
# Run PCA to summarise the major patterns of variation between samples
pca <- prcomp(t(log_cpm))
dim(pca$x)
rownames(pca$x)

# The sample coordinates live in pca$x
# ...so if I ran colnames(pca$x) I get "PC1" "PC2" "PC3" "PC4" "PC5" "PC6"

# Thousands of gene-expression measurements per sample → PCA 
# → a few summary dimensions (PC1, PC2, etc.) describing the major differences between samples

# PC1 → captures the most variation
# PC2 → captures the second most
# PC3 → captures the third most
# PC4 → captures the fourth most

# Look at how much variation each component captures
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# Inspect variation explained by the first four (highest-ranked) PCs
percent_var[1:4]

# Breaking this down...
# prcomp() stores a standard deviation for each component in pca$sdev
# Variance is the square of standard deviation, so ^2
# Dividing each variance by the total gives each component's share
# Multiplying by 100 makes it a percentage

# Create a data frame containing each sample's coordinates
# on PC1 and PC2 (the two components explaining the most variation)
pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# Add sample metadata so PCA points can be identified by condition/stage
pca_df <- merge(pca_df, sample_metadata, by = "sample_id")
pca_df

# ...note, a biologically important treatment effect can still be relatively small and therefore not obvious in those first two PCs

### 10. Building the PCA figure with ggplot2
# The single most informative figure in a routine RNA-seq analysis
# So we create it slowly and thoughtfully

# 1. Data frame, first argument to ggplot()
# 2. Then aesthetic mapping, which column becomes x and y, colour
# ... we'll get an empty panel with correctly scaled axes
ggplot(pca_df, aes(x = PC1, y = PC2))

# 3. Now we give geom argument, what is actually drawn
# Six points
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point()

# 4. Now we map the stage to colour and make points bigger
# Note, if the value comes from a column of your data, it goes inside aes()
# If you picked it yourself, it goes outside
# The same holds for colour, fill, size, shape and alpha
ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 4)

# 5. Set the colours explicitly for reporting
# Name them after the factor levels to the colour always follows the group
# Always use levels() e.g., pca_df$condition to check correct naming of data
ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 4) +
  scale_colour_manual(values = c("a" = "#66c2a5", "b" = "#fc8d62")) +
  theme_bw()

# 10.2 Labels and saving

# Set ggplot object for PCA plot
# paste0() glues strings together without a separator
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 4) +
  scale_colour_manual(values = c("a" = "#66c2a5", "b" = "#fc8d62")) +
  labs(
    title = "PCA of normalised expression",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    colour = "Stage"
  ) +
  theme_bw()

# Draw object
p_pca

# Label the points
# I'm using ggrepel rather than manually setting using vjust
# OPTIONAL ADDITION confirmed by unit chair to be acceptable

p_pca_labels <- p_pca +
  geom_text_repel(
    aes(label = sample_id),
    size = 3,
    show.legend = FALSE
  )

# Display graph
p_pca_labels

# Create dir to save the PNG
dir.create("working/figures", showWarnings = FALSE)

# Save as png using ggsave()
# A resolution of 300 dpi is the usual minimum for a journal figure
ggsave(
  filename = "working/figures/pca_plot_final.png",
  plot = p_pca_labels,
  width = 6,
  height = 5,
  dpi = 300
)

# Say a b sample appeared among the a
# List three possible explanations ranging from a labelling mistake in the metadata to genuine biology
# How would we tell them apart?
# ...my answers are a) metadata labeling mistake in the a ggplot mapping or earlier
# b) one sample didn't ferment as expected
# c) sample might have a technical/sample qual issue, lower seq qual or unusual lib size
# ...contamination, poor alignment, etc. I'd go back to QC

# If PC1 separates your conditions and explains a large share of the variance
# ...the experiment has a strong effect and differential expression will find plenty of genes 
# If the conditions do not separate at all, be cautious about any gene list you produce afterwards
# Look at the data before you test it


#########################################
# Differential expression with DESeq2 ###
#########################################

# Week 9: Differential expression with DESeq2
# Input: filtered raw counts and sample metadata
# Output: a ranked table of differentially expressed genes.
# Comparison: b (45 days) vs a (10 days).

# if cols and rows have drifted apart we must stop and examine
# we set a as the first factor level earlier so it is the reference
stopifnot(all(colnames(counts_filtered) == sample_metadata$sample_id))

### 4. Why raw counts?
# A gene with 10 reads is proportionally far noisier than a gene with 10k reads
# One makes much more of a difference proportionately to ten than ten thousand
# DESeq2 insists on the raw integers and does its own normalisation internally

# Normalise for plotting
# Never normalise before a count-based statistical test
# DESeq2, edgeR and limma-voom all want raw counts and will normalise for you


### 5. Building the DESeqDataSet
# DESeq2 stores the counts, the metadata and the experimental design together in one object
# To ensure the three never drift apart

# 5.1 Build the object

# DESeqDataSetFromMatrix is a class, we pass it arguments
# countData is the gene by sample matrix of raw counts
# colData is the metadata, one row per column of the count matrix, in the same order
# design is a formula... ~ is modeled by, so model each gene's expression as a function of the condition column in the metadata
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_metadata,
  design = ~ condition
)
dds

# 5.2 The design formula matters
# ie for batch experiments you would use design = ~ batch + condition
# The variable of interest goes last in the formula
# Everybg before it is a nuisance var to be adjusted for

# Check what DESeq2 bks the comparison is
# We set a as the first factor level, the reference, in w7
# A positive log fold change will mean higher in the b than the a
# R will set levels alphabetically if we don't explicitly tell it
levels(dds$condition)
resultsNames(dds)

############################
######### Step 1 ###########
############################

### 6. Size factors
# Dividing by the total library size is distorted by a handful of extremely highly expressed genes
# This is is a real risk in yeast where a few genes dominate the transcriptome

# DESeq2 is more clever and instead computes a size factor per sample
# a single number saying how deeply that sample was sequenced compared with a typical sample

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# A size factor near 1 means the sample is close to typical
# Above 1 means it was sequenced more deeply than average
# ... and DESeq2 will scale it down accordingly

# Comparing size factors with raw lib sizes from last week
# Should be correlated but not identical
library_sizes <- colSums(counts_filtered)
round(library_sizes / mean(library_sizes), 3)
round(sizeFactors(dds), 3)

# You can extract the normalised counts if you want to look at them
# Not used for testing
norm_counts <- counts(dds, normalized = TRUE)
head(round(norm_counts, 1))

############################
######### Step 2 ###########
############################

### 7. Dispersion
# Dispersion measures how much a gene's counts vary between replicates of the same condition
# ...beyond what the counting process alone would produce
# It is the parameter that decides whether an observed difference between groups is impressive or ordinary

# With only three replicates per group, estimating a gene's dispersion from its own data is hopeless
# DESeq2 solves this by looking at all genes together
# It fits a curve describing how dispersion typically depends on expression level
# ... then pulls each gene's noisy individual estimate towards that curve
# This is called shrinkage

dds <- estimateDispersions(dds)
plotDispEsts(dds)

# Black points are each gene's own estimate
# The red line is the fitted trend
# Blue points are the final, shrunken values actually used for testing

# Note, points circled in blue are genes with unusually high dispersion
# ...which DESeq2 leaves alone rather than shrinking
# because forcing them down would produce false positives

# The black cloud should scatter around the red line and blue should sit close
# A dispersion plot that does not look like this is a red flag

############################
####### Steps 3 and 4 ######
############################

### 8. Running the test
# The DESeq() function performs size factor estimation
# ...dispersion estimation and model fitting in one call

# We ran the first two steps separately so we could see them
# Running DESeq() now simply completes the job and won't redo work already done

dds <- DESeq(dds)

### DESeq should print our steps:
# using pre-existing size factors
# estimating dispersions
# found already estimated dispersions, replacing these
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

# Always, always check levels before providing arguments
levels(dds$condition)

# Assign and check results
# The second name is the numerator and the third is the reference
# So it is critical we have b before a
res <- results(dds, contrast = c("condition", "b", "a"))
res

# The contrast reads as "condition: b versus a"
# A positive log fold change means higher in the b

############################
##### Reading results ######
############################

### 9. Reading the results table
# Look at the column descriptions that DESeq2 stores with the object

mcols(res)$description

# > mcols(res)$description
# "mean of normalized counts for all samples"      
# "log2 fold change (MLE): condition b vs a"
# "standard error: condition b vs a"        
# "Wald statistic: condition b vs a"        
# "Wald test p-value: condition b vs a"     
# "BH adjusted p-values"      

# Check the Week 9 content for a full breakdown

# Quick overview
summary(res)

### 10. Why the adjusted p-value?

# A p-value threshold of 0.05 means each gene has a 5 per cent chance of looking significant anyway

# Count how many genes clear the raw threshold
# sum(res$pvalue < 0.05, na.rm = TRUE)

# Now the adjusted threshold
# sum(res$padj < 0.05, na.rm = TRUE)

# The na.rm = TRUE argument is needed because DESeq2 sets padj to NA for genes it filtered out as untestable
# Without it, sum() would return NA

# pvalue <0.05
# If this gene were truly unchanged, there is only a 5 per cent chance I would see a difference this big
# Perfectly true for one gene, tested once
# Useless as a rule for choosing a gene list
# ...applying it thousands of times guarantees hundreds of mistakes

# padj < 0.05
# Of all the genes on my final list, about 5 per cent are expected to be false positives
# A statement about the list as a whole
# DESeq2 uses the Benjamini-Hochberg procedure to compute it, controlling the false discovery rate

# Always filter on padj 
# The raw p-value is shown so that you can inspect its distribution, which is a useful diagnostic
# but it is never the basis for a gene list

### 10.1 The p-value histogram
# A histogram of the raw p-values is one of the fastest diagnostics available

res_raw_df <- as.data.frame(res)

ggplot(res_raw_df, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.02, boundary = 0) +
  labs(title = "Distribution of raw p-values", x = "p-value", y = "Number of genes") +
  theme_bw()

# We want a flat carpet aross most of the range
# With a spike near zero
# Flat bit is unchanged genes
# Spike is real signal
# A hist that rises towards 1 or that has bumps in the middle means the model is misspecified, do not trust

### 11. Shrinking the fold changes
# A gene with 4 reads in one group and 1 in the other has a log2 fold change of 2, which looks dramatic and means almost nobg
# Low-count genes produce wildly exaggerated fold changes

# Look at the problem first with an MA plot
# puts expression on the x-axis and fold change on the y-axis
# Every point is a gene

plotMA(res, ylim = c(-6, 6), main = "Before shrinkage")

# Notice the funnel shape on the left
# Weakly expressed genes fan out to extreme fold changes, because they are weakly expressed

# Shrinkage pulls those unreliable estimates towards zero while leaving well-measured genes almost untouched
# when an estimate is noisy, trust it less

# The coef argument needs the name DESeq2 gave the comparison
# Print the names and copy the one we want
resultsNames(dds)

# So we'll use condition_b_vs_a
res_shrunk <- lfcShrink(dds, coef = "condition_b_vs_a", type = "apeglm")
plotMA(res_shrunk, ylim = c(-6, 6), main = "After shrinkage")

# Use the shrunken fold changes for ranking genes, for plotting, and for deciding which genes are biologically interesting
# Use the unshrunken p-values and adjusted p-values for deciding which genes are statistically significant
# Shrinkage changes the effect sizes, not the significance

### 12. Extracting the significant genes

# Turn the shrunken results into a plain data frame
# The third line reorders the columns so that gene_id comes first
res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_shrunk)
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
head(res_df)

# Remove the genes DESeq2 could not test, then sort by adjusted p-value
# !is.na() keeps the rows that have a real adjusted p-value
res_df <- res_df[!is.na(res_df$padj), ]
res_df <- res_df[order(res_df$padj), ]
head(res_df, 10)

# Apply both thresholds
# Statistical significance and a meaningful effect size are different requirements
# A good gene list demands both
# The abs() function takes the absolute value
# so this keeps genes that at least doubled in either direction
# Both conditions must hold

###### ERROR HERE #####
# Should be the line below
# sig_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

# I had to tell R to explicitly run the expression as text
eval(parse(text = "sig_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]"))
nrow(sig_genes)
head(sig_genes)

# Split them by direction
# Remember that positive means higher in the b, because a is the reference
up_in_b <- sig_genes[sig_genes$log2FoldChange > 0, ]
down_in_b <- sig_genes[sig_genes$log2FoldChange < 0, ]
nrow(up_in_b)
nrow(down_in_b)

# How many genes did DESeq2 leave untested
# so how many have padj set to NA?
sum(is.na(res$padj))

# How many genes are significant at < 0.01 
# with an absolute log2 fold change greater than 2?
# This is the same expression as sig_genes above
strict_genes <- res_df[res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 2, ]
nrow(strict_genes)

# A log2 fold change above 2 means at least a fourfold change

### 13. Checking a result with your own eyes

# Take the top gene and plot its normalised counts in each sample
top_gene <- res_df$gene_id[1]
top_gene
plotCounts(dds, gene = top_gene, intgroup = "condition")

# The three a reps should sit cla apart from the three b
# If they overlap, look at metadata closely to rule out labelling error

# Build the same figure in ggplot2 for a publishable version
gene_counts <- plotCounts(dds, gene = top_gene, intgroup = "condition", returnData = TRUE)

# returnData = TRUE argument makes plotCounts() hand back the data frame instead of drawing a base plot
ggplot(gene_counts, aes(x = condition, y = count, colour = condition)) +
  geom_point(size = 4) +
  scale_y_log10() +
  labs(title = top_gene, x = "Stage", y = "Normalised count") +
  theme_bw()

# Repeating ggplot2 figure for the second most significant gene
# Because res_df is sorted by padj, the second row is the second most significant gene
second_gene <- res_df$gene_id[2]
gene_counts2 <- plotCounts(dds, gene = second_gene, intgroup = "condition", returnData = TRUE)

ggplot(gene_counts2, aes(x = condition, y = count, colour = condition)) +
  geom_point(size = 4) +
  scale_y_log10() +
  labs(title = second_gene, x = "Stage", y = "Normalised count") +
  theme_bw()

# And for the top three genes, from the homework
# This is my own code not the unit (AI and my own)
# Using index range from 1 to 3 top genes
top_three_genes <- res_df$gene_id[1:3]
top_three_genes

# Extract normalised counts for the three most significant genes
# lapply(top_three_genes, function(gene) { runs function on each top gene
# gives us a list of three data frames so do.call(rbind, ...) stacks vertically
# it then runs the loop of plotcounts to get norm counts on each gene
top_three_counts <- do.call(
  rbind,
  lapply(top_three_genes, function(gene) {
    
    gene_norm <- plotCounts(
      dds,
      gene = gene,
      intgroup = "condition",
      returnData = TRUE
    )
    
    # Record which gene these counts belong to
    gene_norm$gene <- gene
    
    return(gene_norm)
  })
)

# Set gene factor levels to preserve significance order
top_three_counts$gene <- factor(
  top_three_counts$gene,
  levels = top_three_genes
)

# Plot all three genes
ggplot(top_three_counts,
       aes(x = condition, y = count, colour = condition)) +
  geom_point(size = 4) +
  scale_y_log10() +
  facet_wrap(~ gene) +
  labs(
    x = "Stage",
    y = "Normalised count"
  ) +
  theme_bw()



# Week 10: From a gene list to biology
# Input: DESeq2 results and normalised expression from Weeks 8 and 9.
# Comparison: b (45 days) vs a (10 days).

### 1. Breakdown of script today
# 1. Volcano plot
# Effect size against significance. Which genes changed, and how confidently?

# 2. Heatmap
# Expression patterns across samples. Do the replicates behave alike?

# 3. Annotate
# Turn YDR070C into a name a reader recognises.

# 4. Enrich
# Which biological processes are over-represented in the list?


###########################
### 4. The volcano plot ###
###########################

# A volcano plot shows effect size on the x-axis and stat significance on the y-axis
# Genes that are both strongly changed and confidently detected sit in the top corners

# Sig is plotted as the negative log of the adjusted p-value
# So that small p-values become large numbers
# A very significant gene has a padj near zero

# The logarithm of a number smaller than 1 is negative
# every p-value is smaller than 1
# So log10(0.0001) is -4, not 4
# That's why we need the negative, so the most sig sit at the top of figure
# A padj of 0.1 becomes 1. A padj of 0.01 becomes 2. A padj of 0.0001 becomes 4
# Smaller p-values climb higher up the y-axis

res_df$neg_log10_padj <- -log10(res_df$padj)
head(res_df[, c("gene_id", "log2FoldChange", "padj", "neg_log10_padj")])

# Add a column recording each gene's status
# Same thresholds as last week
# we name the groups after the stages we actually compared

# Each line overwrites the default label for the rows that meet its condition
# Order matters with general ase set first and special aferwards
# Because a is the reference level, a positive fold change means higher at 45 days than at 10 days, 
# We use up in b and down because that is the only claim these labels are entitled to make

res_df$status <- "Not significant"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up in b"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down in b"

# Always check for sensible results
table(res_df$status)

# Plot our scatter plot 
# geom_vline() draws vert lines and geom_hline() horizontal
# Up and down in b is the only real claim these labels are entitled to make
# Because a is the reference level a pos fold change means higher at 45 days than at 10
# names in scale_colour_manual() must match the strings in the status column exactly
# The dashed lines mark the thresholds ie counted as significant

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, colour = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = c(
    "Up in b" = "#d73027",
    "Down in b" = "#4575b4",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  labs(
    title = "Differential expression: b versus a",
    x = "log2 fold change (b / a)",
    y = "-log10 adjusted p-value",
    colour = ""
  ) +
  theme_bw()

# Print plot
p_volcano

### 4.1 Labelling the top genes
# Label only the handful of genes worth naming
# geom_text() is given its own data arg
# A geom can use a different data frame from the one supplied to ggplot()
# this is how you label a subset of points while drawing all of them

top10 <- res_df[order(res_df$padj), ][1:10, ]

# I saved to the var so I could save the png with labels
p_volcano <- p_volcano +
  geom_text(
    data = top10,
    aes(label = gene_id),
    colour = "black",
    size = 2.5,
    vjust = -0.8
  )

# Save it
ggsave("figures/volcano_plot.png", p_volcano, width = 7, height = 6, dpi = 300)

### Challenge 1
# Redraw with stricter 0.01 padj
# How many genes change status?

# Copy the three status lines, change the n
# Give the new column a diff name so we can compare two tables

res_df$status_strict <- "Not significant"
res_df$status_strict[res_df$padj < 0.01 & res_df$log2FoldChange > 1] <- "Up in b"
res_df$status_strict[res_df$padj < 0.01 & res_df$log2FoldChange < -1] <- "Down in b"

# Check
table(res_df$status_strict)
# Compare
table(res_df$status)

# Draw volcano plot with strict genes
p_volcano_strict <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, colour = status_strict)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = c(
    "Up in b" = "#d73027",
    "Down in b" = "#4575b4",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey40") +
  labs(
    title = "Differential expression: b versus a",
    x = "log2 fold change (b / a)",
    y = "-log10 adjusted p-value",
    colour = ""
  ) +
  theme_bw()

# Add helpful gene labels
p_volcano_strict <- p_volcano_strict +
  geom_text(
    data = top10,
    aes(label = gene_id),
    colour = "black",
    size = 2.5,
    vjust = -0.8
  )

# Print plot
p_volcano_strict

# Save strict genes plot
ggsave("figures/volcano_strict_plot.png", p_volcano_strict, width = 7, height = 6, dpi = 300)

###########################
### 5. The heat map #######
###########################

# The volcano plot and the heatmap answer different questions

# Volcanoe plot is one point per gene
# ...samples do not appear at all
# Answers: which genes changed, by how much, and how confidently?

# Heat map
# One row per gene, one column per sample
# Answers: do the replicates wib a group behave alike, and do the two groups look different?

# The volcano plot argues that genes changed
# The heatmap shows that the change is consistent across replicates rather than driven by one odd sample
# A paper usually needs both

# We plot only the significant genes

##########################

# Pull out the significant gene identifiers
# then subset the normalised expression matrix
sig_ids <- res_df$gene_id[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1]
sig_ids <- sig_ids[!is.na(sig_ids)]
length(sig_ids)

# Using the lop_cpm from week 6 or 7 and filtering for sig
heat_matrix <- log_cpm[rownames(log_cpm) %in% sig_ids, ]
# Should have as many rows as sig_ids has entries
dim(heat_matrix)

# This should look like normalised cpm values eg 7.043..
# genes in rows and samples in cols
head(heat_matrix)

# Build an annotation data frame so the heatmap can show which sample belongs to which stage
# The row names must match the column names of the heatmap matrix
annotation_col <- data.frame(Stage = sample_metadata$condition)
rownames(annotation_col) <- sample_metadata$sample_id
annotation_col

### Draw heatmap

# scale = "row" converts each gene's values into z-scores
# so that every gene has a mean of zero and a sd of one across the samples
# shows relative rather than absolute expression
# suppose one gene has values around 15 on the log-cpm scale and another has 2
# without scaling the first gene's row is bright and second's is dark
# neither row shows you anybg about the diff between stages
# with scaling both rows are recentered on their own average
# A heatmap of row-scaled data shows patterns, not expression levels
# A bright red cell means the gene is higher in that sample than in its other samples
# It does not mean the gene is highly expressed
# Consider in the context of the PCA from week 8

# Turn off show_rownames when there are more than about fifty genes

pheatmap(
  heat_matrix,
  scale = "row",
  show_rownames = FALSE,
  annotation_col = annotation_col,
  main = "Differentially expressed genes: b vs a"
)

#####################################
### 6. Annotating the gene names ####
#####################################

# Bioconductor supplies annotation packages that translate common name from systematic
# I.e., FMP16 vs YDR070C

# Check columns of library package
columns(org.Sc.sgd.db)

# Today we'll use common which returns the short name
# and description which returns a sentence describing what the gene does

# the mapIds() function looks up one identifier type and returns another
# keys are the identifiers we have
# keytype says what kind they are
# ORF for yeast systematic names
# column says what we want back
# multiVals = "first" handles the genes that map to more than one name by taking the first
res_df$gene_name <- mapIds(
  org.Sc.sgd.db,
  keys = res_df$gene_id,
  column = "COMMON",
  keytype = "ORF",
  multiVals = "first"
)
head(res_df[, c("gene_id", "gene_name", "log2FoldChange", "padj")], 10)

# Add a description
# some genes will return NA
# many yeast open reading frames simply have no common name
# this is a fact about yeast rather than a code error
# keep the systematic identifier as the primary key
# common names change over time whereas systematic do not
res_df$description <- mapIds(
  org.Sc.sgd.db,
  keys = res_df$gene_id,
  column = "DESCRIPTION",
  keytype = "ORF",
  multiVals = "first"
)
head(res_df$description, 3)

### Challenge 2
# print top 5 sig genes id, common name, and l2fc
# order rows by padj for significance, then choose the columns with a vector
head(res_df[order(res_df$padj), c("gene_id", "gene_name", "log2FoldChange")], 5)

#################################
### 7. Functional enrichment ####
#################################

# Gene Ontology over-representation analysis answers whether we have a random sample of the genome
# or if our listed genes are concentrated in particular biological processes
# For each bio process it asks whether more of our sig genes belong to that process
# ...than would be expected by chance

# Pathway enrichment analysis summarises the long gene list into a shorter
# ...interpret-able list of pathways
# Instead of having a list of a thousand genes you might get 50/60 bio pathways
# You can then check which genes are behind these pathways

# Three ingredients
# 1. Gene list of interest, ie diff expressed genes we want to summarise
# 2. List of background genes, eg all genes in human genome
# 3. List of gene sets, essentially groups of related genes
# ...eg if our list has a lot of genes related to breast cancer
# ...we'd need to tell our code which genes are involved in breast cancer

# We can then compare our list to the background list
# That gives us a list of overrepresented pathways
# e.g., our list of diff expressed genes involves a lot of genes involved in breast cancer

# say we know that ald involves inflammatory processes that involve pro-inflammatory cytokines, e.g., IL-6
# the question would be: is there an association between our genes diff expressed in alcoholic liver disease
# ...versus healthy cells, and IL-6 production
# in other words
# ...is our list of differentially expressed genes enriched with genes involved in IL-6 synthesis pathway

# We build a contingency table
# determine whether the fraction of genes of interest in the pathway is higher 
# ...compared to the fraction of genes outside the pathway, i.e., the background set

### Fishers exact test will give you a p-value for the pathway
# if the p-value is very low you can safely say that your list is overrepresented
# E.g., our gene list is enriched with genes involved in IL-6 production
# IL-6 production is an important pathway in alcoholic liver disease compared to healthy liver cells
# does not mean it's upregulated, we'd need to check gene list
# e.g., are those genes associated with the pos regulation of IL-6 production
# are they positive or negative regulators of that pathway

# multiple testing across pathways mean we will get strong results by chance, ie p value
# we correct this using multiple testing correction
# most commonly used method is Benjamini-Hochberg correction

# the p-value is the probability of seeing at least x n of genes out of n genes in our list
# annotated to a particular gene set term
# given the proportion of background genes that are annotated to that term

# the closer the p-value is to zero
# the more significant the particular gene set term associated with the group of genes is

# A common rule to select DEGs is p-val <0.005 and |FC| >2

### 7.1 The gene universe
# The comparison set is called the universe
# A gene that DESeq2 never tested (filtered out NA) could never have appeared on the list at all
# Counting such genes in the background makes the list look more special than it is
# So the universe is every gene we tested, not every gene in the genome

# Compare the lengths of the gene uni versus sig_ids
# The universe should be thousands of genes and the significant list a few hundred
gene_universe <- res_df$gene_id[!is.na(res_df$padj)]
length(gene_universe)
length(sig_ids)

# The principle does not change with the size of the effect
# State the universe you used and make it the genes you tested 
# A reviewer will ask

### 7.2 Running the analysis

# ont = "BP" restricts the search to biological processes
# pAdjustMethod = "BH" is the same Benjamini-Hochberg correction
# ...because we are testing many hypotheses at once, one per GO term
# GeneRatio is the fraction of your significant genes in that process
# BgRatio is the fraction of the universe in that process
# A term is enriched when the first is much larger than the second

ego <- enrichGO(
  gene = sig_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
head(as.data.frame(ego)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")], 10)

# Visualise top ten terms
# The dot size shows how many genes are in each term
# Colour shows adjusted p-value
# A term with only three genes is rarely worth reporting
dotplot(ego, showCategory = 15) +
  labs(title = "Enriched biological processes: b vs a")

# Enrichment analysis described our gene list
# It shows that our gene list is concentrated in particular annotated processes
# ...more than a random list of the same size would be
# That is a useful summary and a good starting point for a hypothesis
# But its terms are highly redundant

# It does NOT show causation
# A process appearing in the list does not mean it drives the phenotype
# It cannot find a process nobody has annotated
# Poorly studied biology is invisible to this method, by construction

# Treat the output as a summary to be interpreted, not a result to be reported uncritically

### Challenge 3
# Running the analysis on up-regulated genes only

# Create our list
up_ids <- res_df$gene_id[res_df$padj < 0.05 & res_df$log2FoldChange > 1]
up_ids <- up_ids[!is.na(up_ids)]

# Run analysis
# Update the gene arg to point to up_ids
ego_up <- enrichGO(
  gene = up_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
head(as.data.frame(ego_up)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")], 10)

# Visualise top ten terms for upregulated genes
dotplot(ego_up, showCategory = 15) +
  labs(title = "Enriched biological processes: b vs a in upregulated genes")

# And now for down regulated genes
# Create our list
# Change the log2FoldChange to < -1
# Splitting by direction is almost always more informative than one combined list
down_ids <- res_df$gene_id[res_df$padj < 0.05 & res_df$log2FoldChange < -1]
down_ids <- down_ids[!is.na(down_ids)]

# Run analysis
# Update the gene arg to point to down_ids
ego_down <- enrichGO(
  gene = down_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
head(as.data.frame(ego_down)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")], 10)

# Visualise top ten terms for upregulated genes
dotplot(ego_down, showCategory = 15) +
  labs(title = "Enriched biological processes: b vs a in downregulated genes")

### Save our work!
write.csv(res_df, "working/outputs/annotated_results.csv", row.names = FALSE)
write.csv(as.data.frame(ego), "working/outputs/go_enrichment.csv", row.names = FALSE)
saveRDS(heat_matrix, "working/outputs/heat_matrix.rds")

### Session info
sessionInfo()


