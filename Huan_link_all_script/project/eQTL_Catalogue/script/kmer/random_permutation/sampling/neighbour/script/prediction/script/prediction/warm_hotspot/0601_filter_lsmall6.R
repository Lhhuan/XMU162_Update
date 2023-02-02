library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)


setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/warm_hotspot/")
org <- read.table("warm_hotspot_regions_in_pantissue_sorted_merge_win4518.bed.gz",header = F,sep = "\t") %>% as.data.frame()
org$L <- org$V3 - org$V2 
org1 <-filter(org,L>=6)
write.table(org1[,1:3],"warm_hotspot_region_win4518_large_than6.bed",row.names = F, col.names = F,quote =F,sep="\t")
system("gzip warm_hotspot_region_win4518_large_than6.bed")