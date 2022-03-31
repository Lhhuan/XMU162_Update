# library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
# library(Seurat)
# library(clusterProfiler)
# library(org.Hs.eg.db)
library(readr)
# library(EnsDb.Hsapiens.v79)
library(parallel)

setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/prediction/Brain_Cerebellum/output/")
# org <- read_tsv("03_smiliarity.txt",col_names =T)%>%as.data.frame()
org <- fread("03_smiliarity.txt",header =T,sep="\t")%>%as.data.frame()
idfile <-fread("01_merge_all_tissue_cis_sig_eQTL_hotspot_egene_idx.txt",header =T,sep="\t")%>%as.data.frame()
split_point = max(idfile$hotspotidx) + 1 -1

# org1 <- org[which(org>0.098324805)]

# f = mclapply(1:ncol(org), function(i){
#   cat(i,"\n")
#   x = org[,i]
#   ind = which(x>0.098324805)

f = mclapply(1:ncol(org), function(i){
  cat(i,"\n")
  x = org[,i]
  ind = which(x>=0.247)
  if(length(ind)>1){ 
    value = x[ind]
    data.frame(hotspot_id = ind-1, 
              gene_id = i + split_point,
              smiliarity = value,
              stringsAsFactors = FALSE)
  }
}, mc.cores = 1)

aaa <-do.call(rbind,f)

write.table(aaa,"04_transform_and_filter_predict.txt",col.names=T,row.names =F,quote=F,sep="\t")

bbb <-aaa[order(-aaa$smiliarity),]
cutoff <-bbb[round(nrow(bbb)*0.1),"smiliarity"]
top <- filter(bbb,smiliarity >=cutoff)
write.table(top,"04_transform_and_filter_predict_top0.1.txt",col.names=T,row.names =F,quote=F,sep="\t")

