library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(Seurat)
library(R.utils)
library(reshape2)
library(parallel)
library(pheatmap)

setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/prediction/Brain_Cerebellum/output/")
f <- "../data/ENCFF870NPA.h5"
balance_factors <- rhdf5::h5read(f, "balance_factors")
interactions  <- rhdf5::h5read(f, "interactions")
bin_p <- rhdf5::h5read(f, "bin_positions")

p_bin <-t(bin_p)%>%as.data.frame()
chr_bin_range <- rhdf5::h5read(f, "chr_bin_range")
chrs <- rhdf5::h5read(f, "chrs")

chr <-data.frame(id= c(0:24),Chr=chrs)
colnames(p_bin) <-c("id","start","end");
segment <- inner_join(chr,p_bin,by="id")
segment$pos <- paste(segment$Chr,segment$start,segment$end,sep="_")

dfin <- data.frame(interactions)
colnames(dfin) <- segment$pos
rownames(dfin) <- segment$pos
dfin$segment = rownames(dfin)

longin <-melt(dfin,id="segment")
colnames(longin) <- c("segment1","segment2","value")
longin_t <-filter(longin,value !="NaN")
longin_t1 <- filter(longin_t,value >0)

# te <-head(longin_t1,30)
tmp1 <- str_split_fixed(longin_t1$segment1,pattern='_',n=3)
tmp2 <- str_split_fixed(longin_t1$segment2,pattern='_',n=3)
re <-data.frame(chr_1=tmp1[,1],start1=tmp1[,2],end1=tmp1[,3],chr_2=tmp2[,1],start2=tmp2[,2],end2=tmp2[,3],value=longin_t1$value)
re2 <-re[order(re$chr_1,re$start1),]

write.table(re2,"06_Hi_C_result.bed",row.names=F,col.names=F,quote=F,sep="\t")
gzip("06_Hi_C_result.bed")

re3 <-filter(re2,chr_1=="chr1"&chr_2=="chr1")
write.table(re3,"06_Hi_C_result_chr1.bed",row.names=F,col.names=F,quote=F,sep="\t")
gzip("06_Hi_C_result_chr1.bed")


#------------
re <-data.frame(chr_1=tmp1[,1],start1=tmp1[,2],end1=tmp1[,3],chr_2=tmp2[,1],start2=tmp2[,2],end2=tmp2[,3],value=longin_t1$value)
re3 <-filter(re,chr_1=="chr1"&chr_2=="chr1")
re4 <-re3[order(re3$chr_1,re3$start1),]
write.table(re4,"06_Hi_C_result_chr1.bed",row.names=F,col.names=F,quote=F,sep="\t")
gzip("06_Hi_C_result_chr1.bed")

quantile(re$value)
ff <-filter(re,value >2.0281520)
write.table(re4,"06_Hi_C_result_top25.bed",row.names=F,col.names=F,quote=F,sep="\t")
gzip("06_Hi_C_result_chr1.bed")