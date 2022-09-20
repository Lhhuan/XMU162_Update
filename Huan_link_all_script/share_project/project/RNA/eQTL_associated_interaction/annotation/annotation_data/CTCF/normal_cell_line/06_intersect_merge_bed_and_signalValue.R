library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)
setwd("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/CTCF/normal_cell_line/hg38")

GRCh38_chrom_size <- "/share/Projects/huanhuan/ref_data/gencode/GRCh38_chrom.sizes"
chr_set <-paste0("chr",c(1:22))
# for (marker in markers){

peak <-"05_normal_cell_line_ctcf_sort_union_merge.bed.gz"
signalvalue <- "05_normal_cell_line_ctcf_signalvalue_sorted.bed.gz"
merge_signalvalue <- "normal_cell_line_ctcf_merge_peak_interact_signalValue_sorted.bed.gz"
system(paste("bedtools intersect -a",peak,"-b",signalvalue,"-wa -wb |gzip >",merge_signalvalue, sep=" "))
org <- read.table(merge_signalvalue,header = F,sep = "\t") %>% as.data.frame()
colnames(org)[c(1:3,7)] <- c("Chr","start","end","signalvalue")
dat <- org%>%group_by(Chr,start,end)%>%summarize(mean_signalvalue =mean(signalvalue),median_signalvalue =median(signalvalue),max_signalvalue = max(signalvalue))%>%as.data.frame()
fdat <-filter(dat,Chr %in%chr_set)
fmean <-"normal_cell_line_ctcf_merge_mean_signalvalue.bedgraph"
fmedian <-"normal_cell_line_ctcf_merge_median_signalvalue.bedgraph"
fmax <-"normal_cell_line_ctcf_merge_max_signalvalue.bedgraph"
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
