library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)
library(reshape2)
library(R.utils)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/")
Cluster <-read.table("../../11_chr1_louvain_pca10_k300.txt",header = T,sep = "\t") %>% as.data.frame()
org <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/enrichment/hotspot_cutoff_0.176_marker_jaccard_index_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()



chr1_org <- filter(org,hotspot_chr=="chr1")
chr1_org$hotspot <-paste0(chr1_org$hotspot_chr,":",chr1_org$hotspot_start,"-",chr1_org$hotspot_end)
chr1_org<-left_join(chr1_org,Cluster[,c("cluster","hotspot")],by="hotspot")


for(marker in unique(chr1_org$Marker)){
    dat <-filter(chr1_org,Marker==marker)
    p1 <- ggplot(dat, aes(x=as.factor(cluster), y=jaacard_index,group=cluster,fill = as.factor(cluster))) + 
    geom_boxplot(outlier.colour = NA)+
    theme_bw()
    pdf(paste0("16_",marker,"_jaacard_boxplot_cluster.pdf"),height=6,width=7)
    print(p1)
    dev.off()
    print(marker)
}


