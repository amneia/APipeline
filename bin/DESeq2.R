#!/usr/bin/env Rscript
##Expression Quantification
##Rscript DESeq2.R <Counts_Files> <CompareFiles>
## CompareFiles: CompareFiles.csv：Control,Treat

#Input parameter
args <- commandArgs(TRUE)
Counts_Matrix <- args[1]
CompareFile <- args[2]


#Library Packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)


#Input Files
CompareFile <- read.csv(CompareFile, header = TRUE)
Counts_Matrix <- read.csv(Counts_Matrix, header = TRUE)
row.names(Counts_Matrix) <- Counts_Matrix[,1]
Counts_Matrix <- Counts_Matrix[,-1]
Counts_Matrix <- select(Counts_Matrix,append(CompareFile$Control,CompareFile$Treat))


##DEGLists Object Construction
coldata <- data.frame(condition = factor(rep(c('control', 'treat'), each = length(CompareFiles\$Control)), levels = c('control', 'treat')))
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = coldata, design= ~condition)

#FoldChange Calculation
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'treat', 'control'))
res <- data.frame(res,stringsAsFactors = FALSE, check.names = FALSE)

#Differential Expression Gene Filtering
res <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)),]
res[which(res$log2FoldChange >= 1 & res$padj < 0.01),'sig'] <- 'up'
res[which(res$log2FoldChange <= -1 & res$padj < 0.01),'sig'] <- 'down'
res[which(abs(res$log2FoldChange) <= 1 | res$padj >= 0.01),'sig'] <- 'none'
res_select <- subset(res1, sig %in% c('up', 'down'))
write.csv(res1_select, file = "Diff_Genes_DEseq2.csv")
