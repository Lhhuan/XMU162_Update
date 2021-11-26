library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)
types <-c("ALL","chr1_22")
for(type in types){
  type <- "ALL"
  setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/",type))
  org1<-read.table("08_prepare_number_ROC_refine_factor_count.txt",header =  T,sep = "\t") %>% as.data.frame()


  rs<-data.frame()
  for(i in c(1:nrow(org1))){
    data<-c(org1[i,]$Number_of_factor_in_hotspot_TP,
    org1[i,]$Number_of_non_factor_in_hotspot_FP,
    org1[i,]$Number_of_factor_in_non_hotspot_FN,
    org1[i,]$Number_of_non_factor_in_non_hotspot_TN)
    data2<-matrix(data,nrow=2)
    aa <- fisher.test(data2)
    ssss<- aa$p.value <0.05
    tmp<-data.frame(row.names=NULL,Pvalue=round(aa$p.value,3),OR=round(aa$estimate,3),Significant=ssss,Factor=org1[i,]$Factor,Cutoff =org1[i,]$Cutoff)
    # tmp<-data.frame(row.names=NULL,Pvalue=aa$p.value,OR=aa$estimate,Significant=ssss,Type=org1[i,]$Type,Factor=org1[i,]$Factor,Overlap_fraction =org1[i,]$Fraction)
    # tmp=data.frame(Pvalue=aa$p.value,OR=aa$estimate,Type=org1[i,]$Type,Factor=org1[i,]$Factor )
    # rs<-rbind(rs,tmp)
    print(i)
    rs <- bind_rows(rs,tmp)
  }
  write.table(rs,"10_fisher_exact_test_result_factor.txt",row.names = F, col.names = T,quote =F,sep="\t")
}









# intervals = c(6,7,8,9,12,15,18)
# intervals = c(8,12,15,18)
# # interval = 18
# # for(interval in intervals){
# ProcessBedGz <- function(interval = NULL){
#   setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/fisher_exact_test/fisher_result/interval_",interval))
#   org1<-read.table("07_prepare_fisher_test.txt",header =  T,sep = "\t") %>% as.data.frame()


#   rs<-data.frame()
#   for(i in c(1:nrow(org1))){
#     data<-c(org1[i,]$Number_of_factor_in_hotspot,org1[i,]$Number_of_non_factor_in_hotspot,org1[i,]$Number_of_factor_in_non_hotspot,org1[i,]$Number_of_non_factor_in_non_hotspot)
#     data2<-matrix(data,nrow=2)
#     aa <- fisher.test(data2)
#     ssss<- aa$p.value <0.05
#     tmp<-data.frame(row.names=NULL,Pvalue=round(aa$p.value,3),OR=round(aa$estimate,3),Significant=ssss,Type=org1[i,]$Type,Factor=org1[i,]$Factor,Overlap_fraction =org1[i,]$Fraction)
#     # tmp<-data.frame(row.names=NULL,Pvalue=aa$p.value,OR=aa$estimate,Significant=ssss,Type=org1[i,]$Type,Factor=org1[i,]$Factor,Overlap_fraction =org1[i,]$Fraction)
#     # tmp=data.frame(Pvalue=aa$p.value,OR=aa$estimate,Type=org1[i,]$Type,Factor=org1[i,]$Factor )
#     # rs<-rbind(rs,tmp)
#     print(i)
#     rs <- bind_rows(rs,tmp)
#   }
#   write.table(rs,"08_fisher_exact_test_result.txt",row.names = F, col.names = T,quote =F,sep="\t")
# }


# mclapply(intervals,ProcessBedGz, mc.cores = 10)