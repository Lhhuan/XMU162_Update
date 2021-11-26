library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(conflicted)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------
library(Hmisc)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/ROADMAP/Ovary/enhancer/")
    # for (type in types){
        # output_file2<-paste("NHPoisson_emplambda_interval_",j,"_cutoff_",cutoff,"_",tissue,".txt", sep = "")
org<-read.table("01_normal_E097-H3K4me1_H3K27ac.narrowPeaksorted_overlap.gz",header = F,sep = "\t") %>% as.data.frame()
org$length <-org$V3-org$V2
title_name <- "Ovary enhancer"
pdf("Distribution_of_length_of_enhancer.pdf") 
p1<-ggplot(org, aes(x =length)) +geom_histogram(position="identity")+
        xlab("Length") + ylab("Count")+p_theme+ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))

print(p1)
dev.off()




colnames(org)[1] <-"chr"
colnames(org)[2] <-"start"
colnames(org)[3] <-"end"
write.table(org,"02_enhancer_length.txt",row.names = F, col.names = T,quote =F,sep="\t")

mean_length<- round(mean(org$length))
median.length <- median(org$length)
rs <-data.frame(mean_length =mean_length,median_length = median.length)
write.table(rs,"02_enhancer_length_statistics.txt",row.names = F, col.names = T,quote =F,sep="\t")
