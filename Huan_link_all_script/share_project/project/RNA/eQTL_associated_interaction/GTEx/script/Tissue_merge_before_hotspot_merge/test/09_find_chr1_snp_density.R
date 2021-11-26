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

# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/")

org <- read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/tmp_output/08_chr_hotspot_ref_alt_uniformity.bed.gz",header = F,sep = "\t") %>% as.data.frame()

colnames(org) <-c("hotspot_chr","hotspot_start","hotspot_end","eQTL_chr","eQTL_start","eQTL_end","ref","alt")
aaa <-org%>%group_by(hotspot_chr,hotspot_start,hotspot_end)%>%summarize(qtl_count=n())%>%as.data.frame()

aaa$hotspot_length = aaa$hotspot_end -aaa$hotspot_start
aaa$qtl_ratio <-aaa$qtl_count/aaa$hotspot_length
pdf("chr1_qtl_ratio.pdf")
hist(aaa$qtl_ratio)
dev.off()

pdf("chr1_qtl_ratio.pdf")
p1 <- ggplot(aaa,aes(x=1,y=qtl_ratio))+geom_boxplot(width =0.1)+geom_violin()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())
print(p1)
dev.off()

pdf("chr1_qtl_length.pdf")
# hist(aaa$hotspot_length)

p1 <- ggplot(aaa,aes(x=1,y=hotspot_length))+geom_boxplot(width =0.1)+geom_violin()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())
print(p1)
dev.off()
