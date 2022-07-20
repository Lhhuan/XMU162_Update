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

org <- read.table("../../../figure/hit_hospot_ratio_kmer.txt",header = T,sep = "\t") %>% as.data.frame()

kmer_need <-filter(org,ratio>0.1)
setwd("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/sampling/emplambda_0/figure/")
load("all_random_kmer_need_test_value_1_100.Rdata")
all_leng=nrow(all_random_kmer_need_test_value)/100/449
random_kmer_more_than0 <-filter(all_random_kmer_need_test_value,value>0)
random_more0_count <- random_kmer_more_than0 %>% group_by(seq,random_ID)%>%summarise(number=n())%>%data.frame()
random_more0_count$ratio <-random_more0_count$number/all_leng
random <-filter(random_more0_count,seq %in%kmer_need$seq)

colnames(kmer_need)[2:3] <-c("hotspot_count","hotspot_ratio")
dat <-inner_join(random,kmer_need,by="seq")
dat$diff <- dat$ratio -dat$hotspot_ratio 
dat$diff0 <-NA
dat[which(dat$diff>=0),"diff0"]=1
dat[which(dat$diff<0),"diff0"]=0
fdat <-dat%>%group_by(seq)%>%summarise(n=sum(diff0))%>%data.frame()
fdat$p <-fdat$n/100

fdat0 <-filter(fdat,p<0.05)
wilcox <- read.table("hotspot_emplambda_0_random_kmer_wilcox.txt",header = T,sep = "\t") %>% as.data.frame()
wilcox_sig <-filter(wilcox,p_value <0.05)
wilcox_sig_in <-filter(wilcox_sig,seq %in%fdat0$seq)
