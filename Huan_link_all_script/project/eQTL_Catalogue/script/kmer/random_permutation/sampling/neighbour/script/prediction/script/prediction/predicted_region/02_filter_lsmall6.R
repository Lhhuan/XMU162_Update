library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)


setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/predicted_region/")
org <- read.table("predicted_regions_win5000.bed.gz",header = F,sep = "\t") %>% as.data.frame()
org$L <- org$V3 - org$V2 
org1 <-filter(org,L>=6)
write.table(org1[,1:3],"predicted_regions_win5000_large_than6.bed",row.names = F, col.names = F,quote =F,sep="\t")
# system("gzip predicted_regions_win5000_large_than6.bed")