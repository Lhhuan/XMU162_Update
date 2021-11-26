library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(readr)
library(mygene)
setwd("/home/huanhuan/project/kinase/2_strong_gene-impact/script/")
di<-"~/project/kinase/2_strong_gene-impact"


# tfs <- read_tsv(file_path('data/kinase.txt', col_names=F)
tfs<-read.table(file.path(di, "data/DepMap/Achilles_gene_effect.csv"),header = T,sep = ",") %>% as.data.frame()
n<-ncol(tfs)
n