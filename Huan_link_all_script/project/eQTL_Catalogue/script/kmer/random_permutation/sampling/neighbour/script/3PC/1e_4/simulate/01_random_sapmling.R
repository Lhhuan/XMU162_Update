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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/3PC/1e_4/simulate/output/markers/")
org <-read.table("/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()

set.seed(1)
C1 <-org[sample(1:nrow(org), 1000, replace = FALSE),]
org1 <- setdiff(org,C1)
set.seed(1)
C2 <-org1[sample(1:nrow(org1), 1000, replace = FALSE),]

org2 <- setdiff(org1,C2)
set.seed(1)
C3 <-org2[sample(1:nrow(org2), 1000, replace = FALSE),]

write.table(C1,"01_random_c1.bed",col.names=F,row.names=F,quote=F,sep="\t")
write.table(C2,"01_random_c2.bed",col.names=F,row.names=F,quote=F,sep="\t")
write.table(C3,"01_random_c3.bed",col.names=F,row.names=F,quote=F,sep="\t")