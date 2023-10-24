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
group <- rep(c('control', 'treat'), each = (length(CompareFile$Control)))
dgelist <- DGEList(counts = Counts_Matrix, group = group)

#Low Expression Data Filtering
keep <- rowSums(cpm(dgelist) > 1) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

#Normalization
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

#Variance multiplier calculation
design <- model.matrix(~group)

dge <- estimateDisp(dgelist_norm, design, robust = TRUE)

fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
lrt <-  as.data.frame(lrt)

lrt <-  lrt[order(lrt$FDR, lrt$logFC, decreasing = c(FALSE, TRUE)), ]
lrt[which(lrt$logFC >= 1 & lrt$FDR < 0.01),'sig'] <- 'up'

lrt[which(lrt$logFC <= -1 & lrt$FDR < 0.01),'sig'] <- 'down'

lrt[which(abs(lrt$logFC) <= 1 | lrt$FDR >= 0.01),'sig'] <- 'none'

gene_diff_select <- subset(lrt, sig %in% c('up', 'down'))

write.csv(gene_diff_select, file = "Diff_Genes_edgeR.csv")
