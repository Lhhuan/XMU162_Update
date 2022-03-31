
library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(Seurat)
library(R.utils)
library(reshape2)
library(parallel)
library(pheatmap)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/")
org <- read.table("../../output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed",header = F,sep = "\t") %>% as.data.frame()

colnames(org) <-c("chr","start","end")

hot_c <-org%>%group_by(chr)%>%summarise(count=n())%>%as.data.frame()

save(hot_c,file ="hotspot_number_chr.Rdata")

write.table(hot_c,"hotspot_number_chr.txt",row.names=F,quote=F,sep="\t")