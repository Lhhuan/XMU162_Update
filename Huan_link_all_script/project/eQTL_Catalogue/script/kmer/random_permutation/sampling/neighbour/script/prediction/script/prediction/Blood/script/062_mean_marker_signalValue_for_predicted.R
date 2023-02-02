library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)

# markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","CTCF","TFBS","CHROMATIN_Accessibility")
markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","TFBS","CHROMATIN_Accessibility")
# GRCh38_chrom_size <- "/share/Projects/huanhuan/ref_data/gencode/GRCh38_chrom.sizes"
# chr_set <-paste0("chr",c(1:22))
# for (marker in markers){

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/Blood/output/signalValue/mean/")
predcited_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_warm_hotspot_region_win4518_large_than6.bed.gz")
    org <- read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
    colnames(org)[c(1:3,7)] <- c("Chr","start","end","signalvalue")
    dat <- org%>%group_by(Chr,start,end)%>%summarize(mean_signalvalue =mean(signalvalue),median_signalvalue =median(signalvalue),max_signalvalue = max(signalvalue))%>%as.data.frame()
    output_file <-paste0(marker,"_max_mean_median_warm_hotspot_region_win4518_large_than6.bed")
    write.table(dat,output_file,row.names = F, col.names = T,quote =F,sep="\t")
    print(marker)
    system(paste("gzip",output_file,sep=" "))
}

lapply(markers,predcited_signal)
