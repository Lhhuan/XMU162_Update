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
nodes<-read.table("10_nodes_anno_count.txt",header = T,sep = "\t") %>% as.data.frame()
edges<-read.table("10_edges_anno_count.txt",header = T,sep = "\t") %>% as.data.frame()

nodes$markers <-gsub("CHROMATIN_Accessibility","CA",nodes$markers)
pdf("./figure/11_nodes_features.pdf",width = 4,height = 4)
p1 <-ggplot(data = nodes, mapping = aes(x =reorder(markers,fraction) , y = fraction *100)) + geom_bar(stat = 'identity', fill = "#74959A", width=0.6)+
  p_theme+scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black",angle=60,hjust=1,size=8))+
    labs(x="Features",y="Percentage of hotspot(%)")

print(p1)
dev.off()


edges$Types <-gsub("_","-",edges$Types)
edges$Types <-gsub("TSS","PRO",edges$Types)

pdf("./figure/11_edges_features.pdf",width = 2.5,height = 4)
p1 <-ggplot(data = edges, mapping = aes(x =reorder(Types,fraction) , y = fraction *100)) + geom_bar(stat = 'identity', fill = "#74959A", width=0.6)+
  p_theme+scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black",angle=60,hjust=1,size=8))+
    labs(x="Types",y="Percentage of edge(%)")

print(p1)
dev.off()

TAD<-read.table("../output/edges_annotation/09_success_egene_hotspot_TAD.bed.gz",header = F,sep = "\t") %>% as.data.frame()
link_TAD <- TAD[,c(8:11)]%>%unique()
colnames(link_TAD) <-c("CHR","start","end","tissue")
link_TAD$length <-link_TAD$end - link_TAD$start
a <-as.data.frame(table(log10(link_TAD$length)))

a$Var1 <-as.numeric(as.character(a$Var1))
pdf("./figure/11_distribution_of_success_TAD.pdf",width = 4,height = 4)
p1 <-ggplot(data = a, mapping = aes(x = Var1, y = Freq)) + geom_bar(stat = 'identity', fill = "#74959A", width=0.011)+ #
  p_theme+scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))+ #scale_x_continuous(breaks=c(4.5,5,5.5,6,6.5,7))+
    labs(y="Frequency",x="Log10(Length of TAD)")

print(p1)
dev.off()