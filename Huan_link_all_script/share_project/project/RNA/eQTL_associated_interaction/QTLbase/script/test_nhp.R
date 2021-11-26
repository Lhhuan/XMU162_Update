library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(data.table)
library(parallel)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/")
org <- read.table("NHPoisson_emplambda_interval_6cutoff_7.3_all_eQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()


before_org<-read.table("01_merge_QTLbase_all_pop_eQTL_1kg_Completion.txt.gz",header = T,sep = "\t") %>% as.data.frame()



