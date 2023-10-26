#!/usr/bin/env Rscript
##Saturation Curve Plotting
##Rscript Saturation.R <Gene_Number_File> <Read_Number_File>

#Input parameter
library(ggplot2)
Number_File <- commandArgs(TRUE)
Gene_Number_File <- Number_File[1]
Read_Number_File <- Number_File[2]

#Read File
Genes_Number_File <- read.table(Gene_Number_File, sep = "")
colnames(Genes_Number_File) <- c("Sample_Name","Gene_Number")
Reads_Number_File <- read.table(Read_Number_File, sep = "")
colnames(Reads_Number_File) <- c("Sample_Name","Reads_Number")

Count <- cbind(Reads_Number_File,Genes_Number_File$Gene_Number)
colnames(Count)[3] <- "Gene_Number"

#Plotting
a <- max(Count$Reads_Number)
p <- ggplot(Count,aes(x=Reads_Number,y=Gene_Number))+
        geom_line()+
        geom_point()+
        labs(x='Reads number', y='Dectected gene number',title = 'Sequencing saturation analysis')+
        theme(legend.title=element_text(family = "Times New Roman",face = "bold"),
              legend.text=element_text(family = "Times New Roman",face = "bold"),
              panel.background = element_blank(),
              panel.grid.major = element_blank(),
              axis.line = element_line(size = 0.5, colour = 'black'),
              axis.title = element_text(size=15,family = "Times New Roman",face = "bold"),
              axis.text = element_text(size = 12,family = "Times New Roman",face = "bold"),
              plot.title = element_text(size = 20,family = "Times New Roman",face = "bold",hjust = 0.5),
              legend.position = 'right',
              legend.background = element_blank(),
              legend.key = element_blank(),
              legend.key.size = unit(5,'mm'))

ggsave(p,filename = "Saturation.png")
