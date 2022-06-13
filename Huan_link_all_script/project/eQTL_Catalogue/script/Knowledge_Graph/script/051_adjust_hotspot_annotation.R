library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/Knowledge_Graph/output/")
#------------------nodes
nodes_anno <- read.table("./nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt.gz",header = T,sep = "\t") %>% as.data.frame()
nodes_anno <-nodes_anno[,c("Hotspot","H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3","CTCF","CHROMATIN_Accessibility","TFBS","Enhancer")]

for(i in c(2:ncol(nodes_anno))){
    nodes_anno[!is.na(nodes_anno[,i]),i] <-1
    nodes_anno[is.na(nodes_anno[,i]),i] <-0
    print(i)
}

write.table(nodes_anno,"nodes_annotation/hotspot_annotation.txt",row.names = F, col.names = T,quote =F,sep="\t")
