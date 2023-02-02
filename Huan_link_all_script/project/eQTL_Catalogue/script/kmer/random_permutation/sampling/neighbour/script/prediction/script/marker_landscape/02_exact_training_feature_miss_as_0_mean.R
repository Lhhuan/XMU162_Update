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
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/marker_landscape/output/markers_na0/mean/")

# markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","CTCF","TFBS","CA")
markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","TFBS","CA")
# marker="H3K27ac"
read_posi <-function(marker){
    org<-fread(paste0(marker,"_mean_positive.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    rownames(org) <-org[,4]
    org1 <-org[,c(7:106)]
    colnames(org1)<-c(paste(marker,1:100,sep="_"))
    dat <- rowMeans(org1)%>%data.frame()
    colnames(dat) <-marker
    return(dat)
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
    dat <- rowMeans(org1)%>%data.frame()
    colnames(dat) <-marker
    return(dat)
}
tmp_n <- lapply(markers,read_nega)
nega <- do.call(cbind,tmp_n)
nega$hotspot <-rownames(nega)
nega$class <-0

train <-bind_rows(posi,nega)
colnames(train)[9] <- "CHROMATIN_Accessibility"
train <-train[,c("CHROMATIN_Accessibility","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","TFBS","hotspot","class")]

write.table(train,"../../02_training_dataset_miss_as_0_mean.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip ../../02_training_dataset_miss_as_0_mean.txt")

ALL_train <- fread("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/02_whole_genome_training_dataset.txt.gz",header = T,sep = "\t") %>% as.data.frame()
ALL_train1 <-ALL_train[,1:402]
ALL_train2 <- left_join(ALL_train1,train[,1:10],by="hotspot")


lable2<- fread("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/GC_content/gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()

# lable2$cluster <-paste0("C",lable2$cluster)
lable2$hotspot <- paste0(lable2$hotspot_chr,":",lable2$hotspot_start,"-",lable2$hotspot_end)
lable2 <- lable2[,5:6]
colnames(lable2)[1] <- "Lable2"
ALL_train3 <- left_join(ALL_train2,lable2,by="hotspot")
ALL_train3$Lable2[is.na(ALL_train3$Lable2)] <-0
colnames(ALL_train3)[4:5] <-c("segment","Lable1")
ALL_train3 <- ALL_train3[,c(1:5,412,6:411)]

write.table(ALL_train3,"../../02_training_dataset_all_features_marker_landscape_miss_as_0_mean.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip ../../02_training_dataset_all_features_marker_landscape_miss_as_0_mean.txt")

