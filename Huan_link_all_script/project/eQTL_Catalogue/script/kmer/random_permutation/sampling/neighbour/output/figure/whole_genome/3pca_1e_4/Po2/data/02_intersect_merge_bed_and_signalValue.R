library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/Po2/data/")

# markers = c("Human_FACTOR","HISTONE_MARK_AND_VARIANT","Human_CHROMATIN_Accessibility")
GRCh38_chrom_size <- "/share/Projects/huanhuan/ref_data/gencode/GRCh38_chrom.sizes"
chr_set <-paste0("chr",c(1:22))
# for (marker in markers){
# signal<-function(marker= NULL){
    marker="PolII_hg38_chr1_22"
    org <- read.table("PolII_hg38_chr1_22_sorted_merge_bed_signal.bed.gz",header = F,sep = "\t") %>% as.data.frame()
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
    print(marker)
# }
# marker="/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/Human_FACTOR/POL2"
# lapply(marker,signal)