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
# for(j in c(10:17)){
ProcessBedGz <- function(j = NULL){
    file_name<-paste("../../output/ALL_eQTL/NHPoisson_emplambda_interval_",j,"cutoff_7.3_all_eQTL.txt.gz",sep="")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org2<-filter(org,emplambda != "NA")

    # emplambda<-org2$emplambda
    result <-ks.test(org2$emplambda,"pgamma", 3, 2)
    tmp<-data.frame(Pvalue=result$p.value,D = result$statistic,j=j)
    return(tmp)
    gc()
}



# file2<-files[1:2]
# interval <-c(10:17)
interval <-c(18,20,24,28,32)
test_result <-mclapply(interval, ProcessBedGz, mc.cores = 10)

rs <-data.frame()
for (i in c(1:5)){
    rs<-bind_rows(rs,test_result[[i]])
}
colnames(rs)[3]<-"interval"
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/")
write.table(rs,"ks_test_result.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(aa,"ks_test_result.txt",row.names = F, col.names = T,quote =F,sep="\t")