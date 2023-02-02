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

org <- read.table("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/normal_tissue/output/02_tissue_level_Marker_source.txt",header = T,sep = "\t") %>% as.data.frame()
org <- org[,2:3]%>%unique()
org <-filter(org,marker!="HISTONE_MARK_AND_VARIANT")
histone_markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3")
histone_info <- filter(org,marker %in% histone_markers)
ca_info <- filter(org, !(marker %in% histone_markers))
#==================
GRCh38_chrom_size <- "/share/Projects/huanhuan/ref_data/gencode/GRCh38_chrom.sizes"

chr_set <-paste0("chr",c(1:22))
# for (marker in markers){
for (i in 1:nrow(histone_info)){
    tissue <-histone_info[i,"refine_tissue_label2"]
    marker <- histone_info[i,"marker"]
    path <-paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/",tissue,"/marker/histone_hg38/")
    setwd(path)
    peak <- paste0(marker,"_merge_pos_info_narrow_peak_sorted_merge.bed.gz")
    signalvalue <- paste0(marker,"_merge_pos_info_narrow_peak_signalValue_sorted.bed.gz")
    merge_signalvalue <- paste0(marker,"_merge_peak_interact_signalValue_sorted.bed.gz")
    system(paste("bedtools intersect -a",peak,"-b",signalvalue,"-wa -wb |gzip >",merge_signalvalue, sep=" "))
    org <- read.table(merge_signalvalue,header = F,sep = "\t") %>% as.data.frame()
    colnames(org)[c(1:3,7)] <- c("Chr","start","end","signalvalue")
    dat <- org%>%group_by(Chr,start,end)%>%summarize(mean_signalvalue =mean(signalvalue),median_signalvalue =median(signalvalue),max_signalvalue = max(signalvalue))%>%as.data.frame()
    fdat <-filter(dat,Chr %in%chr_set)
    fmean <-paste0(marker,"_merge_mean_signalvalue.bedgraph")
    fmedian <-paste0(marker,"_merge_median_signalvalue.bedgraph")
    fmax <-paste0(marker,"_merge_max_signalvalue.bedgraph")
    write.table(fdat[,c(1:3,4)],fmean,row.names = F, col.names = F,quote =F,sep="\t")
    write.table(fdat[,c(1:3,5)],fmedian,row.names = F, col.names = F,quote =F,sep="\t")
    write.table(fdat[,c(1:3,6)],fmax,row.names = F, col.names = F,quote =F,sep="\t")
    #==================TO bw
    fmean_bw <-gsub("bedgraph","bw",fmean)
    fmedian_bw <-gsub("bedgraph","bw",fmedian)
    fmax_bw <-gsub("bedgraph","bw",fmax)
    system(paste("bedGraphToBigWig",fmean,GRCh38_chrom_size,fmean_bw,sep =" "))
    system(paste("bedGraphToBigWig",fmedian,GRCh38_chrom_size,fmedian_bw,sep =" "))
    system(paste("bedGraphToBigWig",fmax,GRCh38_chrom_size,fmax_bw,sep =" "))

    gzip(fmean)
    gzip(fmedian)
    gzip(fmax)
    print(c(tissue,marker))
}

for (i in 1:nrow(ca_info)){
    tissue <-ca_info[i,"refine_tissue_label2"]
    marker <- ca_info[i,"marker"]
    path <-paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/",tissue,"/marker")
    setwd(path)
    peak <- paste0(marker,"_merge_pos_info_narrow_peak_sorted_merge.bed.gz")
    signalvalue <- paste0(marker,"_merge_pos_info_narrow_peak_signalValue_sorted.bed.gz")
    merge_signalvalue <- paste0(marker,"_merge_peak_interact_signalValue_sorted.bed.gz")
    system(paste("bedtools intersect -a",peak,"-b",signalvalue,"-wa -wb |gzip >",merge_signalvalue, sep=" "))
    org <- read.table(merge_signalvalue,header = F,sep = "\t") %>% as.data.frame()
    colnames(org)[c(1:3,7)] <- c("Chr","start","end","signalvalue")
    dat <- org%>%group_by(Chr,start,end)%>%summarize(mean_signalvalue =mean(signalvalue),median_signalvalue =median(signalvalue),max_signalvalue = max(signalvalue))%>%as.data.frame()
    fdat <-filter(dat,Chr %in%chr_set)
    fmean <-paste0(marker,"_merge_mean_signalvalue.bedgraph")
    fmedian <-paste0(marker,"_merge_median_signalvalue.bedgraph")
    fmax <-paste0(marker,"_merge_max_signalvalue.bedgraph")
    write.table(fdat[,c(1:3,4)],fmean,row.names = F, col.names = F,quote =F,sep="\t")
    write.table(fdat[,c(1:3,5)],fmedian,row.names = F, col.names = F,quote =F,sep="\t")
    write.table(fdat[,c(1:3,6)],fmax,row.names = F, col.names = F,quote =F,sep="\t")
    #==================TO bw
    fmean_bw <-gsub("bedgraph","bw",fmean)
    fmedian_bw <-gsub("bedgraph","bw",fmedian)
    fmax_bw <-gsub("bedgraph","bw",fmax)
    system(paste("bedGraphToBigWig",fmean,GRCh38_chrom_size,fmean_bw,sep =" "))
    system(paste("bedGraphToBigWig",fmedian,GRCh38_chrom_size,fmedian_bw,sep =" "))
    system(paste("bedGraphToBigWig",fmax,GRCh38_chrom_size,fmax_bw,sep =" "))

    gzip(fmean)
    gzip(fmedian)
    gzip(fmax)
    print(c(tissue,marker))
}
