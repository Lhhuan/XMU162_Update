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

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/")
org<-read.csv("6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(org)[1] <-"hotspot"
org_hp <- read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_0.176.bed.gz",header=F,sep="\t")%>%as.data.frame()

org_hp$hotspot <-  paste0(">",org_hp$V1,":",org_hp$V2,"-",org_hp$V3)

kmer_h <-org$hotspot  %>%as.data.frame()
colnames(kmer_h)[1] <-"hotspot"
kmer_h$count <- "a"

diff_h1 <-anti_join(kmer_h,org_hp,by= "hotspot")
diff_h <-anti_join(org_hp,kmer_h,by= "hotspot")