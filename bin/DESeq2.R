#!/usr/bin/env Rscript
##Expression Quantification
##Rscript Gene_Counts.R <Counts_Files>

#Input parameter
args <- commandArgs(TRUE)
Counts_Matrix <- args[1]
CompareFile <- args[2]

#Library Packages
