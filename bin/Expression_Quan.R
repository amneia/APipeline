#!/usr/bin/env Rscript
##Expression Quantification
##Rscript Gene_Counts.R <Counts_Files>

#Input parameter
Counts_Files <- commandArgs(TRUE)

#Results Construction
Counts_Matrix <- as.data.frame(matrix(nrow=0, ncol=2))
colnames(Counts_Matrix) <- c("Sample_Name", "Gene_Number")
for ( i in 1:length(Counts_Files) ){
        Sample_Name <- substr(strsplit(Counts_Files[i], "[.]")[[1]][1],1,(nchar(strsplit(Counts_Files[i], "[.]")[[1]][1])-7))
        Counts_File <- read.table(Counts_Files[i], sep = "\t", header = T)
        Gene_Count <- 0
        for ( j in 1:length(Counts_File[,7]) ){
                if ( Counts_File[j,7] != 0 ){
                       Gene_Count <- Gene_Count + 1
                }
        }
        Count_Matrix <- data.frame(Sample_Name, Gene_Count)
        Counts_Matrix <- rbind(Counts_Matrix,Count_Matrix)
}

write.csv(Counts_Matrix, file = "Gene_Counts.csv",row.names = FALSE )
