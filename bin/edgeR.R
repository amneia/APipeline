#!/usr/bin/env Rscript
##Expression Quantification
##Rscript edgeR.R <Counts_Files> <CompareFiles>
## CompareFiles: CompareFiles.csv：Control,Treat

#Input parameter
args <- commandArgs(TRUE)
Counts_Matrix <- args[1]
CompareFile <- args[2]

#Library packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
library(dplyr)
library(edgeR)

