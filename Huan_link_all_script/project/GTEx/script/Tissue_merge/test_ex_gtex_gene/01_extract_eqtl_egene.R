
# source("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/OmicCircos_circos_refine_revise_color.R")
library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(tidyverse)
library(reshape2)
library(Seurat)
library(circlize)
library(tibble)


setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/test_ex_gtex_gene/")


org <- read.table("/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL_hg19/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",header = T,sep = "\t") %>% as.data.frame()
#-----
f_org <-filter(org,pval_nominal<5E-8)
f_n_org <-f_org%>%select(Chr,Pos,gene_id)%>%unique()
f_n_org  <- add_column(f_n_org, end=f_n_org$Pos+1, .after = 2)
f_n_org$gene_id <-stringr::str_split_fixed(f_n_org$gene_id,"\\.",2)[,1]
f_n_org$Chr <-paste0("chr",f_n_org$Chr)


write.table(f_n_org,"R_Whole_Blood_cis_sig_eQTL_egene.txt",col.names=F,row.names=F,quote=F,sep="\t")