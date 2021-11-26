library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)


org<-read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/merge_QTL_all_QTLtype_pop.txt.gz",header = T,sep = "\t") %>% as.data.frame()
#----mkdir and set pwd
cutoff = 7.3
#----------------------------------------------
xQTL ="eQTL"
org_pop <- filter(org,QTL_type == xQTL)
i =1
org2<-filter(org_pop, SNP_chr==i)
org2$SNP_pos<-as.numeric(as.character(org2$SNP_pos))
org2$SNP_chr<-as.numeric(as.character(org2$SNP_chr))
#-------------------------------------- snp has unique p，min_p
org_p <- org2%>%dplyr::select(SNP_pos,Pvalue)%>%unique()
# rm(org2)
org_pg <- group_by(org_p, SNP_pos)
org1<-summarise(org_pg,min_p = min(Pvalue))%>%as.data.frame()
# rm(org_pg)
#-----------------------------------
org1 <- org1[order(org1$SNP_pos), ]
org1$chr <-"1"
org1$start<-org1$SNP_pos
org1$end <- org1$start +1
org1$value <- 1
org1<-org1[,-c(1:2)]
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/")
write.table(org1,"QTLbase_chr1_eQTL.bed",row.names = F, col.names = F,quote =F,sep="\t")
print("11")