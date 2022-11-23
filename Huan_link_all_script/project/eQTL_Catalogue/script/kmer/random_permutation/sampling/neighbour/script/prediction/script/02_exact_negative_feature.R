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
library(mclust)
library(umap)
library(uwot)
library(kernlab)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
org<-read.csv("01_all_negative_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure")
# load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/emplambda_0/06_permutation_test_1000_sig_kmer.Rdata")
load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/06_permutation_test_1000_sig_kmer.Rdata")
sigK <-fdat0

colnames(org)[1] <-"hotspot"
org$hotspot <-gsub(">","",org$hotspot)
rownames(org)=org$hotspot
org <-org[,-1]

Sorg <-org[,which(colnames(org) %in% sigK$seq)]
