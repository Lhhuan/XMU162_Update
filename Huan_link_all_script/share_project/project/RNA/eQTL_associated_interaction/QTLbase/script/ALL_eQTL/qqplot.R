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

ProcessBedGz <- function(j = NULL){
     file_name<-paste("../../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz",sep="")
     org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
     org2<-filter(org,emplambda != "NA")

     emplambda<-org2$emplambda
     title =paste("Q-Q Plot for ","gamma"," distribution: interval ",j,sep ="" )
     figure_name1<-paste("./add_qqplot/Q_Q_Plot_for_gamma_distribution_interval_",j,"_add0",".png",sep ="")
     png(figure_name1)
     emplambda[which(emplambda==0)]=emplambda[which(emplambda==0)]+1e-12
     qualityTools::qqPlot(emplambda,"gamma",main=title)
     dev.off()
     print("add")
     #-----------------------------------------
     figure_name2<-paste("./remove_qqplot/Q_Q_Plot_for_gamma_distribution_interval_",j,"_remove0",".png",sep ="")
     png(figure_name2)
     emplambda<-org2$emplambda
     # emplambda[which(emplambda==0)]=emplambda[which(emplambda==0)]+1e-12
     emplambda<-emplambda[which(emplambda!=0)]
     qualityTools::qqPlot(emplambda,"gamma",main=title)
     dev.off()
     print("remove")
     print(j)
}



# # file2<-files[1:2]
# interval <-c(10:17)
# # interval <-c(18,20,24,28,32)
# test_result <-mclapply(interval, ProcessBedGz, mc.cores = 10)

interval <-c(12,18,20,24,28,32)
# interval <-c(10:17)
mclapply(interval, ProcessBedGz, mc.cores = 10)