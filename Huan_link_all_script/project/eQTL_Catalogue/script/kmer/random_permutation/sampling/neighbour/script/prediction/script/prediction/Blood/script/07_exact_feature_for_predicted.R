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

nega<-fread("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/warm_hotspot/07_predicted_warm_hotspot_region_gc_kmerdataset.txt.gz",header = T,sep = "\t") %>% as.data.frame()
#========================================
markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","TFBS","CHROMATIN_Accessibility")
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/Blood/output/signalValue/mean/")
nega_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_max_mean_median_warm_hotspot_region_win4518_large_than6.bed.gz")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org <-org[,c(1:4)]
    org$marker <-marker
    return(org)
}

tmp2 <-lapply(markers,nega_signal)
rs2<-do.call(rbind,tmp2)
nega_S <-reshape2::dcast(rs2[,c("Chr","start","end","marker","mean_signalvalue")], Chr+start+end~marker)
nega_S$hotspot <-paste0(nega_S$Chr,":",nega_S$start,"-",nega_S$end)
# nega_S[is.na(nega_S)] <-0
nega <-left_join(nega,nega_S[,c(4:13)],by="hotspot")
nega[is.na(nega)] <-0
# nega <-nega%>%dplyr::select(-Class)
#========================================

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/Blood/output/")
write.table(nega,"07_predicted_warm_hotspot_region_dataset.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip 07_predicted_warm_hotspot_region_dataset.txt")
