library(TCGAbiolinks)
library(SummarizedExperiment)
library(pheatmap)
library(tidyverse)
library(maftools)
library(MCPcounter)
library(devtools)
install.packages("ggplot2")
library(ggplot2)

GDC <- getGDCprojects()
query_TCGA <- GDCquery(project = 'TCGA-CHOL',
                       experimental.strategy = 'RNA-Seq',
                       data.category = 'Transcriptome Profiling',
                       access = "open",
                       sample.type = "Primary Tumor")
output <- getResults(query_TCGA)
GDCdownload(query_TCGA)
tcga_chol_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
chol_matrix <- assay(tcga_chol_data, 'fpkm_unstrand')


# Load probesets and genes data
probesets <- read.table("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt", sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
genes <- read.table("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)

# Filter the 111 genes
filter_gene <- (genes[,4])
#clean up the dot and everything after it
rownames(chol_matrix) <- sub("\\..*", "", rownames(chol_matrix))
colnames(chol_matrix) <- sub("01A-11R-A41I-07", "", colnames(chol_matrix))

# Filter out rows 
filtered_matrix <- chol_matrix[rownames(chol_matrix) %in% filter_gene, ]

min_non_zero <- min(filtered_matrix[filtered_matrix > 0.0000])
filtered_matrix[filtered_matrix == 0.0000] <- min_non_zero

#log2
log2_matrix <- log2(filtered_matrix)
MCPScore=MCPcounter.estimate(log2_matrix,featuresType="ENSEMBL_ID")

#rearrange

desired_order <- rev(c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils", "Endothelial cells", "Fibroblasts"))

MCPScore <- MCPScore[desired_order,]

heatmap((MCPScore),col=colorRampPalette(c("blue","white","red"))(100),Rowv=NA, Colv = TRUE)

tMCPScore <- t(MCPScore)
# Perform PCA
pca_result <- prcomp(tMCPScore)
# Create a PCA plot with PC1 on X-axis and PC2 on Y-axis
plot(pca_result$x[, 1], pca_result$x[, 2], 
     xlab = "PC1", ylab = "PC2", 
     main = "PCA Plot: PC1 vs PC2", 
     pch = 16, col = "blue")
                     