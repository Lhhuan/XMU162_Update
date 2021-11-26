library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(circlize)
library(CMplot)
library(reshape2)
library(tidyverse)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")
org<-fread("../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL_allQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()
aa<- org %>% group_by(chr)
max_pos<-aa %>% summarise(max(t))%>%as.data.frame
min_pos<-aa %>% summarise(min(t))%>%as.data.frame
ss<-cbind(min_pos,max_pos)
write.table(ss,"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/find_max_min_emplambda.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(max_pos,"/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/find_max_emplambda.txt",row.names = F, col.names = F,quote =F,sep="\t")