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
load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/06_permutation_test_1000_sig_kmer.Rdata")
sigK <-fdat0

colnames(org)[1] <-"hotspot"
org$hotspot <-gsub(">","",org$hotspot)
rownames(org)=org$hotspot
# org <-org[,-1]

Norg <-org[,which(colnames(org) %in% sigK$seq)]
Norg$hotspot <- rownames(Norg)
nega_gc <- read.csv("01_all_negative.gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(nega_gc) <- c("Chr","start","end","GC_content")
nega_gc$hotspot <-paste0(nega_gc$Chr,":",nega_gc$start,"-",nega_gc$end)
nega_gc$Class <-0
nega_gc <- nega_gc[,c(1:3,5,6,4)]
nega <-left_join(nega_gc,Norg,by="hotspot")
# nega$Class <-0

#=============positive
load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
Sorg$hotspot <- rownames(Sorg)
posi_gc <- read.csv("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/GC_content/gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(posi_gc)[1:4] <-c("Chr","start","end","GC_content")
posi_gc$hotspot <-paste0(posi_gc$Chr,":",posi_gc$start,"-",posi_gc$end)
posi_gc$Class <-1
posi_gc <- posi_gc[,c(1:3,6,7,4)]

posi <-left_join(posi_gc,Sorg,by="hotspot")

ALL_train <- rbind(nega,posi)
write.table(ALL_train,"02_whole_genome_training_dataset.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip 02_whole_genome_training_dataset.txt")
Tchr1 <- filter(ALL_train,Chr=="chr1")
write.table(Tchr1,"02_chr1_training_dataset.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip 02_chr1_training_dataset.txt")
