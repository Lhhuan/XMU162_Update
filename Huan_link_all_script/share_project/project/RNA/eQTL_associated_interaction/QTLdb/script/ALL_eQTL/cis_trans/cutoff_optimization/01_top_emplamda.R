library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/NHP/cis_1MB")

org<-read.table("NHPoisson_emplambda_interval_18_cutoff_7.3_cis_1MB_eQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org_value<-dplyr::filter(org,emplambda !="NA")
org_value_order <-org_value[order(-org_value$emplambda),]
all_n <-nrow(org)
aa <-seq(0.01,0.1,0.01)
rs<-data.frame()
for(a in aa){
    cutoff <-round(all_n *a)
    value <-org_value_order[cutoff,1]
    print(value)
    b<-nrow(dplyr::filter(org_value_order,emplambda>=value))
    ratio = b/all_n
    tmp<-data.frame(top_fraction_expect=a,emplambda = value,top_number = b,top_fraction_fact=ratio)
    rs<-bind_rows(rs,tmp)
}

write.table(rs,"top_emplambda.txt",c = T, quote =F,sep="\t")
top <-