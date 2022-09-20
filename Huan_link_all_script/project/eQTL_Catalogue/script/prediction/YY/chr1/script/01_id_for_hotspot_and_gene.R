library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)


setwd("/home/huanhuan/project/eQTL_Catalogue/script/prediction/YY/chr1/output/")
org <- read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org) <-c("h_chr","h_start","h_end","gene")
dat <-filter(org,h_chr=="chr1")
dat$hotspot <-paste(dat$h_chr,dat$h_start,dat$h_end,sep="_")
nodes <-unique(c(dat$gene,dat$hotspot))

set.seed(1)
nodes <-nodes[sample(1:length(nodes),length(nodes),replace = FALSE)]
x <-length(nodes)-1
idx <-data.frame(name=nodes,idx =0:x)

write.table(idx,"01_hotspot_egene_idx.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip 01_hotspot_egene_idx.txt")