library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(readr)
library(mygene)
library(R.utils)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/DGLlinker/train/")



# tfs <- read_tsv(file_path('data/kinase.txt', col_names=F)
org<-read.table("IntAct v2020-11-06.csv",header = T,sep = ",") %>% as.data.frame()
set.seed(123)
org1 <-org[sample(nrow(org),nrow(org)),]
org1$Target.node.name <-org$Source.node.name
write.table(org1,"adjust_gene_IntAct.csv",quote = F,sep = ",",row.names = F)
