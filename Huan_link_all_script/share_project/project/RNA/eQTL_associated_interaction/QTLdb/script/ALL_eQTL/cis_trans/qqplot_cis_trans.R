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
    title =paste("Q-Q Plot for ","gamma"," distribution: interval ",j," ",name, sep ="" )
    figure_name1<-paste("./",name,"/add_qqplot/Q_Q_Plot_for_gamma_distribution_interval_",j,"_add0",".png",sep ="")
    png(figure_name1)
    emplambda[which(emplambda==0)]=emplambda[which(emplambda==0)]+1e-12
    qualityTools::qqPlot(emplambda,"gamma",main=title)
    dev.off()
    print("add")
    #-----------------------------------------
    figure_name2<-paste("./",name,"/remove_qqplot/Q_Q_Plot_for_gamma_distribution_interval_",j ," ",name,"_remove0",".png",sep ="")
    png(figure_name2)
    emplambda<-org2$emplambda
    # emplambda[which(emplambda==0)]=emplambda[which(emplambda==0)]+1e-12
    emplambda<-emplambda[which(emplambda!=0)]
    qualityTools::qqPlot(emplambda,"gamma",main=title)
    dev.off()
    print("remove")
    print(j)
}




names <-c("cis_10MB","cis_1MB","trans_10MB","trans_1MB")
interval <-c(10:17)
mclapply(names, function(name){
    print(name)
    mclapply(interval, function(num){
        ProcessBedGz(j = num, name = name)
    }, mc.cores = 8)
#---------------
}, mc.cores = 5) 