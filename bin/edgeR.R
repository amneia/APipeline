#!/usr/bin/env Rscript
##Differential Expression Analysis(edgeR)
##Rscript edgeR.R <Counts_Files> <CompareFiles>

#Input parameter
args <- commandArgs(TRUE)
Counts_Matrix <- args[1]
CompareFile <- args[2]

#Library packages
library(dplyr)
