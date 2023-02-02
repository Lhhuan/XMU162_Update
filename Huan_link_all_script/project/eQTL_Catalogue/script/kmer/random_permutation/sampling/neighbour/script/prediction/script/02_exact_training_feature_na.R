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
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
org<-fread("01_all_negative_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
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


#========================================
markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","CTCF","TFBS","CHROMATIN_Accessibility")
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/signalValue/mean/")
nega_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_max_mean_median_1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz")
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
nega <-left_join(nega,nega_S[,c(4:14)],by="hotspot")
# nega[is.na(nega)] <-0
#========================================



#=============positive
load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
Sorg$hotspot <- rownames(Sorg)
posi_gc <- read.csv("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/GC_content/gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(posi_gc)[1:4] <-c("Chr","start","end","GC_content")
posi_gc$hotspot <-paste0(posi_gc$Chr,":",posi_gc$start,"-",posi_gc$end)
posi_gc$Class <-1
posi_gc <- posi_gc[,c(1:3,6,7,4)]

posi <-left_join(posi_gc,Sorg,by="hotspot")

#===========================
setwd("/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/signalValue/hotspot/mean")
hotspot_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_max_mean_median_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org <-org[,c(1:4)]
    org$marker <-marker
    return(org)
}
tmp2 <-lapply(markers,hotspot_signal)
rs2<-do.call(rbind,tmp2)
posi_S <-reshape2::dcast(rs2[,c("Chr","start","end","marker","mean_signalvalue")], Chr+start+end~marker)
posi_S$hotspot <-paste0(posi_S$Chr,":",posi_S$start,"-",posi_S$end)
posi <-left_join(posi,posi_S[,c(4:14)],by="hotspot")
# posi[is.na(posi)] <-0
#==========================
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
ALL_train <- rbind(nega,posi)
write.table(ALL_train,"02_whole_genome_training_dataset_na.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip 02_whole_genome_training_dataset_na.txt")
Tchr1 <- filter(ALL_train,Chr=="chr1")
write.table(Tchr1,"02_chr1_training_dataset_na.txt",row.names = F, col.names = T,quote =F,sep="\t")
system("gzip 02_chr1_training_dataset_na.txt")
