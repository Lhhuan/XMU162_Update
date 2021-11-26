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
# library(fitdistrplus)
library(qualityTools)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_eQTL/cis_trans")
# for(j in c(10:17)){
ProcessBedGz <- function(j = NULL, name=NULL){
    file_name<-paste("../../../output/ALL_eQTL/cis_trans/NHP/",name,"/NHPoisson_emplambda_interval_",j,"_cutoff_7.3_",name,"_eQTL.txt.gz",sep="")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org2<-filter(org,emplambda != "NA")

    emplambda<-org2$emplambda
    result <-ks.test(emplambda,"pgamma", 3, 2)
    print(j)
    tmp<-data.frame(Pvalue=result$p.value,D = result$statistic,j=j,data_type = name)
    return(tmp)
    gc()
}




names <-c("cis_10MB","cis_1MB","trans_10MB","trans_1MB")
interval <-c(10:17)
mclapply(names, function(name){
    print(name)
    test_result <- mclapply(interval, function(num){
        ProcessBedGz(j = num, name = name)
    }, mc.cores = 8)
#---------------
}, mc.cores = 2) 

