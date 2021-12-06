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

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/")

tissue_merge <- read.table("../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(tissue_merge)<-c("chr","start","end")
org_chr1 <- filter(tissue_merge,chr=="chr1")
org_chr1$value <-1
org_chr1_density <-genomicDensity(org_chr1, window.size=1e4)
chr1_high <-org_chr1_density[which(org_chr1_density$pct==1),]
all_density <-genomicDensity(tissue_merge, window.size=1e4)
all_high <-all_density[which(all_density$pct==1),]

write.table(chr1_high,"../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/05_chr1_density_high.txt",row.names=F,quote=F,sep="\t")
write.table(all_high,"../../output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/05_all_density_high.txt",row.names=F,quote=F,sep="\t")