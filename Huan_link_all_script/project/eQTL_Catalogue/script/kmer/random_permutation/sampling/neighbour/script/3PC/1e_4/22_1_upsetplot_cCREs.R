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
library(ggsci)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/Cis_Regulatory_Elements/")
org <-read.table("hotspot_cluster_cCREs.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)<-c("CHR","start","end","cluster","cre_chr","cre_start","cre_end","name","score","strand","thickStart","thickEnd","reserved","ccre","encodeLabel","zScore","ucscLabel","accessionLabel","description")
org$hotspot <-paste(org$CHR,org$start,org$end,sep=":")
# org$cluster <-as.factor(org$cluster)
org <- org[,c("hotspot","cluster","ucscLabel","name")]
org1 <- unique(org[,c("hotspot","cluster","ucscLabel")])
org1$value <-1
all_hotspot <- read.table("../GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(all_hotspot) <-c("CHR","start","end","cluster")
all_hotspot$hotspot <- paste(all_hotspot$CHR,all_hotspot$start,all_hotspot$end,sep=":")
all_hotspot <- all_hotspot[,c("hotspot","cluster")]

dat <- left_join(all_hotspot,org1,by=c("cluster","hotspot"))

fdat <-dcast(dat,hotspot+cluster~ucscLabel)
fdat[is.na(fdat)] <-0

library(UpSetR)
pdf("22_1_upset_plot_cCREs_all_hotspot.pdf",width=13,height=7)
upset(fdat,mainbar.y.label="Number of hotspot",sets.x.label="Number of hotspot",sets = c("CTCF","enhD","enhP","K4m3","prom"))
dev.off()

f4 <-filter(fdat,cluster==4)
pdf("22_1_upset_plot_cCREs_cluster4.pdf",width=13,height=7)
upset(f4,mainbar.y.label="Number of hotspot",sets.x.label="Number of hotspot",sets = c("CTCF","enhD","enhP","K4m3","prom"))
dev.off()
