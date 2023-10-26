#!/usr/bin/env Rscript
##Percentage of each sample reads after decoding
##Rscript Decode_Results.R <Fastq_Multx_Result_1> <Fastq_Multx_Result_2>

#Input parameter
Fastq_Multx_Result <- commandArgs(TRUE)
Fastq_Multx_Result_1 <- Fastq_Multx_Result[1]
Fastq_Multx_Result_2 <- Fastq_Multx_Result[2]

#Input Files
Fastq_Multx_Result_1 <- read.table(Fastq_Multx_Result_1, sep = "")
Fastq_Multx_Result_2 <- read.table(Fastq_Multx_Result_2, sep = "")

#Calculation of the proportion
Decode_Pro <- paste(sprintf("%.2f",(Fastq_Multx_Result_1$Id + Fastq_Multx_Result_2$Id)*100/(sum(Fastq_Multx_Result_2$Id)*2)),"%",sep="")
Sample_Name <- row.names(Fastq_Multx_Result_2)

#Output Results
Decode_Results <- data.frame(Sample_Name, Decode_Pro)
write.csv(Decode_Results,file = "./Decode_Results.csv")
