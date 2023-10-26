#!/usr/bin/env Rscript
##Saturation Curve Plotting(All)
##Rscript Saturation.R <Gene_Number_File> <Read_Number_File>

#Input parameter
library(ggplot2)
Number_Files <- commandArgs(TRUE)

#Read File
Gene_Number_Files <- subset(Number_Files,grepl("Gene_Number",Number_Files))
Reads_Number_Files <- subset(Number_Files,grepl("Reads_Number",Number_Files))

Genes_Number <- read.table(Gene_Number_Files[1], sep = "")
colnames(Genes_Number) <- c("Sample_Name","Splited_Name","Gene_Number")
Reads_Number <- read.table(Reads_Number_Files[1], sep = "")
colnames(Reads_Number) <- c("Splited_Name","Reads_Number")

for ( i in 2:length(Gene_Number_Files)){
  Genes_Number_File <- read.table(Gene_Number_Files[i], sep = "")
  colnames(Genes_Number_File) <- c("Sample_Name","Splited_Name","Gene_Number")
  Genes_Number <- rbind(Genes_Number, Genes_Number_File)
  Reads_Number_File <- read.table(Reads_Number_Files[i], sep = "")
  colnames(Reads_Number_File) <- c("Splited_Name","Reads_Number")
  Reads_Number <- rbind(Reads_Number, Reads_Number_File)
}

Counts <- merge(Genes_Number,Reads_Number,by = "Splited_Name")

#Plotting
a <- max(Counts$Reads_Number)
p <- ggplot(Counts,aes(x=Reads_Number,y=Gene_Number,color=Sample_Name))+
        geom_line(aes(group=Sample_Name))+
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
