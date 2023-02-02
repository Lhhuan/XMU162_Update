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
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/cancer/script/")
#=========================kmer
all_kmer <- fread("../../predicted_region/03_predicted_regions_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/06_permutation_test_1000_sig_kmer.Rdata")
sigK <-fdat0

colnames(all_kmer)[1] <-"hotspot"
all_kmer$hotspot <-gsub(">","",all_kmer$hotspot)
rownames(all_kmer)=all_kmer$hotspot
kmer <-all_kmer[,which(colnames(all_kmer) %in% sigK$seq)]
kmer$hotspot <- rownames(kmer)
#============================gc 
gc <- read.csv("../../predicted_region/03_predicted_regions_gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(gc) <- c("Chr","start","end","GC_content")
gc$hotspot <-paste0(gc$Chr,":",gc$start,"-",gc$end)
#============================kmer gc
gc_kmer <- left_join(gc[,c("hotspot","GC_content")],kmer,by="hotspot")
tmp_hotspot <- gc_kmer[,c("hotspot","GC_content")]
#============================
ca_tf <- read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/082_final_ca_tf_TCGA_sample.txt",header = T,sep = "\t") %>% as.data.frame()
histone <- read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/092_final_histone_TCGA_mark.txt",header = T,sep = "\t") %>% as.data.frame()

org <- bind_rows(ca_tf,histone)
org <-unique(org[,1:4])
org$tmp_tissue <- paste(org$TCGA,org$cell_line_refine,sep="/")


# org <- org[c(1:3,8:10),]
histone_markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3")

mean_bin <- function(marker=NULL){
    marker_signal <- fread(paste0(pre_dir,"/",marker,"_mean.bed.gz"),header = F,sep = "\t")%>% as.data.frame()
    rownames(marker_signal) <-marker_signal[,4]
    marker_signal1 <-marker_signal[,c(7:106)]
    # colnames(marker_signal1)<-c(paste(marker,1:100,sep="_"))
    dat <- rowMeans(marker_signal1)%>%data.frame()
    marker_refine <- gsub("Human_CHROMATIN_Accessibility","CHROMATIN_Accessibility", marker)
    marker_refine <- gsub("Human_FACTOR","TFBS", marker_refine)
    colnames(dat) <-marker_refine
    #=========================有些marker 某些片段是NA
    dat$hotspot <- rownames(dat)
    dat1 <-left_join(tmp_hotspot,dat,by="hotspot")
    dat2 <- data.frame(dat1[,3])
    rownames(dat2) <- dat1$hotspot
    colnames(dat2) <- marker_refine
    return(dat2)
}

all_marker<-c("CHROMATIN_Accessibility","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3","TFBS")
length(unique(org$tmp_tissue))
# for (i in 61:120){
for (i in 92:120){
    tissue =unique(org$tmp_tissue)[i]
    print(paste0(tissue," start"))
    pre_dir <- paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/cancer/cell_line/",tissue,"/marker/markers_na0/mean")
    T_info <- filter(org,tmp_tissue==tissue)
    tmp_n <- lapply(T_info$factor,mean_bin)
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
    fname <- paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/cancer/cell_line/",tissue,"/marker_matrix_0_mean_allFeatute.txt")
    write.table(fdat,fname,row.names = F, col.names = T,quote =F,sep="\t")
    gzip(fname)
    # print(tissue)
    print(paste0(tissue," end"))
}

