# Load necessary libraries
library(DESeq2)
library(dplyr)
library(ggplot2)

# Set working directory
setwd("/05152025_rrnaseq_DPM/")

# Read count data
data <- read.delim("combine1.txt", header = TRUE, check.names = FALSE)
data <- na.omit(data)  # Remove any rows with NA

# Prepare count matrix
data2 <- data[,-1]                 # Drop gene name column for count matrix
rownames(data2) <- data[,1]       # Use gene names as rownames

# Filter genes: keep those expressed in at least 4 samples
keep <- rowSums(data2 > 0) >= 4
filtered_counts <- data2[keep,]
filtered_counts <- as.matrix(filtered_counts)

# Create sample metadata
coldata <- data.frame(
  sample = colnames(data2),
  condition = factor(rep(c("B", "A"), each = 3))  # Adjust as needed
)
rownames(coldata) <- coldata$sample

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = coldata,
  design = ~ condition
)

# Run DESeq2
dds <- DESeq(dds)

# Export results table
res <- results(dds)
write.table(res, file = "DPM2VSD12.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Export normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.table(norm_counts, file = "normalized_counts.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Plot expression of a gene (example: "SON")
plotCounts(dds, gene = "SON", returnData = TRUE) %>%
  ggplot(aes(x = condition, y = count)) +
  geom_bar(stat = 'summary', fill = 'grey', color = "black", position = "dodge", width = 0.5) +
  geom_jitter(fill = "white", pch = 21, size = 3, stroke = 1, alpha = 1) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 25, colour = 'black', angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20, colour = 'black'),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )
