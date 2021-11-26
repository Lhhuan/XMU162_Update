library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)
# types <-c("ALL","chr1_22")

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_select_exclude_eQTL/Cis_eQTL/ROC/interval_18/ALL/split_non_factor/")
org1<-read.table("08_prepare_number_ROC_factor_count.txt",header =  T,sep = "\t") %>% as.data.frame()


rs<-data.frame()
for(i in c(1:nrow(org1))){
  data<-c(org1[i,]$Number_of_factor_in_hotspot_TP,
  org1[i,]$Number_of_non_factor_in_hotspot_FP,
  org1[i,]$Number_of_factor_in_non_hotspot_FN,
  org1[i,]$Number_of_non_factor_in_non_hotspot_TN)
  data2<-matrix(data,nrow=2)
  aa <- fisher.test(data2)
  ssss<- aa$p.value <0.05
  tmp<-data.frame(row.names=NULL, Tissue = org1[i,]$Tissue, Random_number = org1[i,]$Random_number, Pvalue=round(aa$p.value,3),OR=round(aa$estimate,3),Significant=ssss,Factor=org1[i,]$Factor,Cutoff =org1[i,]$Cutoff)
  # tmp<-data.frame(row.names=NULL,Pvalue=aa$p.value,OR=aa$estimate,Significant=ssss,Type=org1[i,]$Type,Factor=org1[i,]$Factor,Overlap_fraction =org1[i,]$Fraction)
  # tmp=data.frame(Pvalue=aa$p.value,OR=aa$estimate,Type=org1[i,]$Type,Factor=org1[i,]$Factor )
  # rs<-rbind(rs,tmp)
  print(i)
  rs <- bind_rows(rs,tmp)
}
write.table(rs,"10_fisher_exact_test_result_factor.txt",row.names = F, col.names = T,quote =F,sep="\t")





