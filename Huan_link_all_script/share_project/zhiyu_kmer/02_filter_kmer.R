library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)
library(reshape2)
library(parallel)
library(tibble)

setwd("/home/huanhuan/zhiyu_kmer/")

org<-read.csv("6mers_uc_us_no_log.csv",header = T,sep = ",") %>% as.data.frame()
colnames(org)[1] <-"hotspot"
hotspot <-read.table("~/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/figure/hotspot_original_random_kmer_wilcox.txt",header = T,sep = "\t") %>% as.data.frame() 
sig_hotspot <-filter(hotspot,p_value<0.05)
power <-org[,colnames(org) %in% sig_hotspot$seq]
power_f <-add_column(power, nodes=org$hotspot, .before = 1)
write.table(power_f,"hotspot_sig_kmer_value.txt",row.names = F, col.names = T,quote =F,sep="\t")