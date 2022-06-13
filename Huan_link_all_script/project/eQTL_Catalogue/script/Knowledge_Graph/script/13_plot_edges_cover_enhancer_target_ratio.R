library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(tidyverse)
library(reshape2)
library(Seurat)
library(circlize)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
                                                
setwd("/home/huanhuan/project/eQTL_Catalogue/script/Knowledge_Graph/output/")
org<-read.table("./edges_annotation/12_cover_enhancer_target_type_ratio.txt.gz",header = T,sep = "\t") %>% as.data.frame()

org$type <-gsub("enhancer_target","Enhancer_Target",org$type)
pdf("./figure/13_edges_features_enhancer_type.pdf",width = 3,height = 4)
p1 <-ggplot(data = org, mapping = aes(x =reorder(type,ratio) , y = ratio *100)) + geom_bar(stat = 'identity', fill = "#74959A", width=0.6)+
  p_theme+scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black",angle=60,hjust=1,size=8))+
    labs(x="Types",y="Percentage of interactions(%)")

print(p1)
dev.off()
