library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(parallel)
library(conflicted)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/NHP/cis_10MB/")


interval7<-read.table("NHPoisson_emplambda_interval_7_cutoff_7.3_cis_10MB_eQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()

org2<-dplyr::filter(interval7, chr=="19")
max(org2$t)
#63762728




all_cis10MB <- read.table("../../cis_10MB_eQTL_pop_1kg_Completion.txt.gz",header =T,sep = "\t")%>% as.data.frame()
