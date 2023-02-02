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
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/marker_landscape/output/markers_na/mean/")

markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","TFBS","CA")
# marker="H3K27ac"
read_posi <-function(marker){
    org<-fread(paste0(marker,"_mean_positive.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    rownames(org) <-org[,4]
    org1 <-org[,c(7:106)]
    colnames(org1)<-c(paste(marker,1:100,sep="_"))
    return(org1)
}

tmp_p <- lapply(markers,read_posi)
posi <- do.call(cbind,tmp_p)
posi$hotspot <-rownames(posi)
posi$class <-1

read_nega <-function(marker){
    org<-fread(paste0(marker,"_mean_negative.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    rownames(org) <-org[,4]
    org1 <-org[,c(7:106)]
    colnames(org1)<-c(paste(marker,1:100,sep="_"))
    return(org1)
}
tmp_n <- lapply(markers,read_nega)
nega <- do.call(cbind,tmp_n)
nega$hotspot <-rownames(nega)
nega$class <-0

train <-bind_rows(posi,nega)

write.table(train,"../../02_training_dataset_miss_as_na.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip ../../02_training_dataset_miss_as_na.txt")
