library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)
library(data.table)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/normal_tissue/")
#=========================kmer
all_kmer <- fread("../predicted_region/03_predicted_regions_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/06_permutation_test_1000_sig_kmer.Rdata")
sigK <-fdat0

colnames(all_kmer)[1] <-"hotspot"
all_kmer$hotspot <-gsub(">","",all_kmer$hotspot)
rownames(all_kmer)=all_kmer$hotspot
kmer <-all_kmer[,which(colnames(all_kmer) %in% sigK$seq)]
kmer$hotspot <- rownames(kmer)
#============================gc 
gc <- read.csv("../predicted_region/03_predicted_regions_gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(gc) <- c("Chr","start","end","GC_content")
gc$hotspot <-paste0(gc$Chr,":",gc$start,"-",gc$end)
#============================kmer gc
gc_kmer <- left_join(gc[,c("hotspot","GC_content")],kmer,by="hotspot")
#============================
org <- read.table("./output/02_tissue_level_Marker_source.txt",header = T,sep = "\t") %>% as.data.frame()
org <- org[,2:3]%>%unique()
org <-filter(org,marker!="HISTONE_MARK_AND_VARIANT")
# org <- org[c(1:3,8:10),]
histone_markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3")

mean_bin <- function(marker=NULL){
    if(marker%in%histone_markers){
        # fname <-paste0(pre_dir,"/histone_hg38/markers_na0/mean/",marker,"_",level,"bed.gz")
        marker_signal <- fread(paste0(pre_dir,"/histone_hg38/markers_na0/mean/",marker,"_mean.bed.gz"),header = F,sep = "\t")%>% as.data.frame()
        rownames(marker_signal) <-marker_signal[,4]
        marker_signal1 <-marker_signal[,c(7:106)]
        # colnames(marker_signal1)<-c(paste(marker,1:100,sep="_"))
        dat <- rowMeans(marker_signal1)%>%data.frame()
        colnames(dat) <-marker
        return(dat)
    }else{
        marker_signal <- fread(paste0(pre_dir,"/markers_na0/mean/",marker,"_mean.bed.gz"),header = F,sep = "\t")%>% as.data.frame()
        rownames(marker_signal) <-marker_signal[,4]
        marker_signal1 <-marker_signal[,c(7:106)]
        # colnames(marker_signal1)<-c(paste(marker,1:100,sep="_"))
        dat <- rowMeans(marker_signal1)%>%data.frame()
        marker_refine <- gsub("Human_CHROMATIN_Accessibility","CHROMATIN_Accessibility", marker)
        marker_refine <- gsub("Human_FACTOR","TFBS", marker_refine)
        colnames(dat) <-marker_refine
        return(dat)
    }
}
all_marker<-c("CHROMATIN_Accessibility","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","TFBS")

for (i in 1:length(unique(org$refine_tissue_label2))){
    tissue=unique(org$refine_tissue_label2)[i]
    print(paste0(tissue," start"))
    pre_dir <- paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/",tissue,"/marker")
    T_info <- filter(org,refine_tissue_label2==tissue)
    tmp_n <- lapply(T_info$marker,mean_bin)
    marker_mean <- do.call(cbind,tmp_n)
    n=length(all_marker)
    if(ncol(marker_mean)<n){
        miss_markers <- setdiff(all_marker,colnames(marker_mean))
        m=ncol(marker_mean)+1
        marker_mean[,m:n] <-NA 
        colnames(marker_mean)[m:n] <-miss_markers
    }
    marker_mean$hotspot <-rownames(marker_mean)
    marker_mean <-marker_mean[,c("CHROMATIN_Accessibility","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","TFBS","hotspot")]
    fdat <- left_join(gc_kmer,marker_mean,by="hotspot")
    fname <- paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/",tissue,"/marker_matrix_0_mean_allFeatute.txt")
    write.table(fdat,fname,row.names = F, col.names = T,quote =F,sep="\t")
    gzip(fname)
    # print(tissue)
    print(paste0(tissue," end"))
}

