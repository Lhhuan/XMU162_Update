library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

pwd_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script"
setwd(pwd_dir)
org<-read.table("../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz",header = T,sep = "\t") %>% as.data.frame()
n_all <-nrow(org) #293623477
cis_1M<-filter(org,cis_trans_1MB=="cis")
n_cis_1M <-nrow(cis_1M) #271386276
n_trans_1M =n_all-n_cis_1M #22237201

cis_10M<-filter(org,cis_trans_10MB=="cis")
n_cis_10M <-nrow(cis_10M) # 275994248
n_trans_10M =n_all-n_cis_10M #17629229