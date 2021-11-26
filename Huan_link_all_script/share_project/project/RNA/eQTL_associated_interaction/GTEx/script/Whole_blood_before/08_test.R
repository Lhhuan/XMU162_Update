library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/NHP/")

org2<-read.table("NHPoisson_emplambda_interval_18_cutoff_7.3_Whole_Blood.txt.gz",header = T,sep = "\t") %>% as.data.frame()

aa <-org2[which(org2$emplambda==1)]
bb <-dplyr::filter(org2,emplambda==1)
chr1 <-dplyr::filter(org2,chr==1)
gz1 <- gzfile("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/NHP/interval_18_chr1.txt.gz", "w")
write.table(chr1,gz1,row.names = F, col.names = T,quote =F,sep="\t")
close(gz1)

interval = 18


