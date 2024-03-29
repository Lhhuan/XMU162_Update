library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/hotspot_associate_study_and_resource/output/")
org <- read.csv("adjust_tissue_need_study_for_hotspot.tsv",header = T,sep = "\t") %>% as.data.frame()
anno <- read.csv("adjust_tissue_need_study_for_hotspot_annotation.tsv",header = T,sep = "\t") %>% as.data.frame()
refine_anno <-left_join(org[,1:2],anno,by=c("study","qtl_group"))
write.table(refine_anno,"01_adjust_tissue_need_study_for_hotspot_annotation.tsv",row.names = F, col.names = T,quote =F,sep="\t")
