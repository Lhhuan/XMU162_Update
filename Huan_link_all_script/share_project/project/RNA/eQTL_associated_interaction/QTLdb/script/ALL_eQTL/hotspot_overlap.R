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
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_eQTL/")



# j=6

# file_name<-paste("../../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz",sep="")
# org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()


ProcessBedGz <- function(j = NULL){
    file_name<-paste("../../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz",sep="")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org2<-filter(org,emplambda != "NA")
    rm(org)
    org2$key<-paste(org2$chr,org2$t,sep=":")
    org2<-org2[,-c(2:3)]
    colnames(org2)[1] <- j
    if (j==19){
        hotspot <-filter(org2,emplambda >0.4)
        return(hotspot)
    }else{
        hotspot <-filter(org2,emplambda >0.2)
        return(hotspot)
    }
    rm(org2)
}




# interval <-c(12,18,20,24,28,32)

interval <-c(6,7,8,9,12,15,19)
read_file <- mclapply(interval, ProcessBedGz, mc.cores = 7)

tmp<-read_file[[1]]

for (i in length(read_file)){
    tmp<-left_join(tmp,read_file[[i]]) 
}

