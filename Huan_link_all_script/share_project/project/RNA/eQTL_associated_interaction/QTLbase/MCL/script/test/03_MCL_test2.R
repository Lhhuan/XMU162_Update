#------------compute-0-6 /state/partition1/huan

library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/MCL/output")

LD<-read.table("/state/partition1/huan/02_hotspot_in_1kg_phase3_v5_eur_1MB_10000_0.2.ld.gz",header = T,sep = "") %>% as.data.frame()

difference<-LD$CHR_A - LD$CHR_B
unique_difference <-unique(difference)
#0
LD<-LD[,-c(1,2,4,5)]
#---------------------count snp
SNP_A <-unique(LD$SNP_A)
SNP_A<-data.frame(snp = SNP_A)
SNP_B <-unique(LD$SNP_B)
SNP_B<-data.frame(snp = SNP_B)
all_snp <-rbind(SNP_A,SNP_B)
n<-nrow(all_snp) #17188253
#17188253 #10000

#-----------------------make matrix for MCL
LD_M <-dcast(LD,SNP_A~SNP_B)
rownames(LD_M)<-LD_M$Var1
LD_M<-LD_M[,-1]
LD_M[is.na(LD_M)] <- 0
#---------------------start MCL
library(MCL)
result <-mcl(mat, addLoops = T, expansion = 2, inflation = 1.4)

#-------------






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