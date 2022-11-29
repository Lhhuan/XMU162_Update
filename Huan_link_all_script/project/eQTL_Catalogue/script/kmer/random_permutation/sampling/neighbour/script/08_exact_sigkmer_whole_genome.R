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
library(data.table)
setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/")
org<-fread("6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure")
# load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/emplambda_0/06_permutation_test_1000_sig_kmer.Rdata")
load("../06_permutation_test_1000_sig_kmer.Rdata")
sigK <-fdat0

colnames(org)[1] <-"hotspot"
org$hotspot <-gsub(">","",org$hotspot)
rownames(org)=org$hotspot
org <-org[,-1]

Sorg <-org[,which(colnames(org) %in% sigK$seq)]
save(Sorg,file ="08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
# save(Sorg,file ="08_permutation_wilocx_overlap_sig_kmer_0_1000_V3.4.Rdata",version=3.4)
write.table(Sorg,"08_permutation_wilocx_overlap_sig_kmer_0_1000.txt",row.names = F, col.names = T,quote =F,sep="\t")
library(kernlab)
pca1<- kpca(~.,data=Sorg, kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)