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

interval <-c(6,7,8,9,12,15,19)
aa<-data.frame()
for(j in interval){
# ProcessBedGz <- function(j = NULL){
    file_name<-paste("../../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz",sep="")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org2<-filter(org,emplambda != "NA")

    # emplambda<-org2$emplambda
    # aa<-data.frame()
    for(alpha in c(1:3)){
        for(beta in c(1:3)){
            result <-ks.test(org2$emplambda,"pgamma", alpha, beta)
            tmp<-data.frame(Pvalue=result$p.value,D = result$statistic,interval=j,alpha= alpha,beta =beta)
            aa<-bind_rows(aa,tmp)
            # gc()
        }
    }
    # return(aa)
}

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/")
write.table(aa,"ks_test_result.txt",row.names = F, col.names = T,quote =F,sep="\t")
print("aaa")
#---------------------------


interval <-c(6,7,8,9,12,15,19)
test_result <-mclapply(interval, ProcessBedGz, mc.cores = 10)

rs <-data.frame()
for (i in length(test_result)){
    rs<-bind_rows(rs,test_result[[i]])
}
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/")
write.table(rs,"ks_test_result.txt",row.names = F, col.names = T,quote =F,sep="\t")
print("aaa")