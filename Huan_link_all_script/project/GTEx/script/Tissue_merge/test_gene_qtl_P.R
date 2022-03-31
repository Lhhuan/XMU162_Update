library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

# setwd("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/enrichment/figure/0_0.176/")
# org<-read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/01_merge_all_tissue_cis_eQTL_gene.txt.gz",header = T,sep = "\t") %>% as.data.frame()

# thyroid <- read.table("/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL/Thyroid.v8.signif_variant_gene_pairs.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# Nerve <- read.table("/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL/Nerve_Tibial.v8.signif_variant_gene_pairs.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# skin <- read.table("/share/data0/GTEx/data/GTEx_Analysis_v8_eQTL/Skin_Sun_Exposed_Lower_leg.v8.signif_variant_gene_pairs.txt.gz",header = T,sep = "\t") %>% as.data.frame()


# tmp <-bind_rows(thyroid,Nerve,skin)
org<-read.table("/share/data0/GTEx/data/all_tissue_eQTL_hg19.txt.gz",header = T,sep = "\t") %>% as.data.frame()

aa <- filter(org, chr==1 & Pos=="173488286")
aa <- filter(org, chr==1 & Pos=="173492282")

173492282