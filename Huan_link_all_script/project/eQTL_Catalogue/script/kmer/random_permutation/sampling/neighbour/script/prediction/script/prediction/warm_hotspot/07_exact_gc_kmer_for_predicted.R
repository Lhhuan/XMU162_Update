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
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/warm_hotspot/")
org<-fread("06_predicted_warm_hotspot_regions_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/06_permutation_test_1000_sig_kmer.Rdata")
sigK <-fdat0

colnames(org)[1] <-"hotspot"
org$hotspot <-gsub(">","",org$hotspot)
rownames(org)=org$hotspot
# org <-org[,-1]

Norg <-org[,which(colnames(org) %in% sigK$seq)]
Norg$hotspot <- rownames(Norg)
nega_gc <- read.csv("06_predicted_warm_hotspot_regions.gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(nega_gc) <- c("Chr","start","end","GC_content")
nega_gc$hotspot <-paste0(nega_gc$Chr,":",nega_gc$start,"-",nega_gc$end)
nega_gc <- nega_gc[,c(1:3,5,4)]
nega <-left_join(nega_gc,Norg,by="hotspot")

# setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
write.table(nega,"07_predicted_warm_hotspot_region_gc_kmerdataset.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip 07_predicted_warm_hotspot_region_gc_kmerdataset.txt")




#========================================
markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","CTCF","TFBS","CHROMATIN_Accessibility")
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/signalValue/mean/")
nega_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_max_mean_median_warm_region_win4518_large_than6.bed.gz")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org <-org[,c(1:4)]
    org$marker <-marker
    return(org)
}

tmp2 <-lapply(markers,nega_signal)
rs2<-do.call(rbind,tmp2)
nega_S <-dcast(rs2[,c("Chr","start","end","marker","mean_signalvalue")], Chr+start+end~marker)
nega_S$hotspot <-paste0(nega_S$Chr,":",nega_S$start,"-",nega_S$end)
# nega_S[is.na(nega_S)] <-0
nega <-left_join(nega,nega_S[,c(4:14)],by="hotspot")
nega[is.na(nega)] <-0
nega <-nega%>%dplyr::select(-Class)
#========================================

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
write.table(nega,"07_predicted_warm_region_dataset.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip 07_predicted_warm_region_dataset.txt")
