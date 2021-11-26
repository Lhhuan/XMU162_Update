library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output")

all_QTL<-fread("./ALL_QTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_QTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()

i=1
aa<-filter(all_QTL,chr==i)
pos<-aa$t
r <- dist(pos,diag = TRUE,method =  "manhattan",upper =T)







qtls = c("eQTL","sQTL","eQTL","mQTL","pQTL","hQTL","caQTL")
for(qtl in qtls){
    filename<-paste("./ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",qtl,".txt.gz",sep="")
    qtl<-fread("./ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_eQTL.txt.gz",header = F,sep = "\t") %>% as.data.frame()
}

ALL_QTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_QTL.txt.gz
eQTL<-fread("./ALL_eQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_eQTL.txt.gz",header = F,sep = "\t") %>% as.data.frame()
mQTL<-fread("./ALL_mQTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_mQTL.txt.gz",header = F,sep = "\t") %>% as.data.frame()



#-----------------------矩阵转换
sss<-dcast(cc,Var1~Var2)
rownames(sss)<-sss$Var1
sss<-sss[,-1]
sss[is.na(sss)] <- 0