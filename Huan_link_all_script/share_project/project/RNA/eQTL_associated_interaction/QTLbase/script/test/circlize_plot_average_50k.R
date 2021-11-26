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
org<-fread("../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_50kb_bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(org)[1] <- "chr"
colnames(org)[2] <- "start1"
colnames(org)[3] <- "end1"
colnames(org)[5] <- "start2"
colnames(org)[6] <- "end2"
colnames(org)[7] <- "xQTL"
colnames(org)[8] <- "emplambda"
org_need<-org[,c('chr','start1','end1','start2','end2','xQTL','emplambda')]
org_need$region<-paste(org_need$start1,org_need$end1,sep="_")
org_a<-org_need%>%group_by(xQTL,chr,region)%>%summarise(mean(emplambda))%>%as.data.frame()
colnames(org_a)[4] <- "mean_emplambda"
data<- spread(org_a,xQTL,mean_emplambda)
nrow(data)
#53583
