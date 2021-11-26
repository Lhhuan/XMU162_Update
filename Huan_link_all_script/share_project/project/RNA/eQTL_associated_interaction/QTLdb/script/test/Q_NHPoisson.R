library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

# R >= 3.1.1 is needed

if (!require("NHPoisson")) {
    install.packages("NHPoisson")
    library("NHPoisson")
}

# setting a random seed for reproducibility
set.seed(1138)



setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")

org<-read.table("../output/chr22_eur.txt",header = T,sep = "\t") %>% as.data.frame()
# dateB <- org$SNP_pos
# BarEv <- POTevents.fun(T = -log10(org$Pvalue), thres = 7, date = dateB)
# Px <- BarEv$Px
# inddat <- BarEv$inddat

# png("emplambda_all_qtl_chr22.png",height = 800,width = 1000)
# emplambdaB <- emplambda.fun(posE = Px, inddat = inddat, t = org$SNP_pos, lint = 500, tit = "EUR_all_chr22")
# dev.off()
i=22
indx <-which(org$ QTL_type == "eQTL")
org1 <- org[indx, ][1:10000, ]
org1 <- org1[order(org1$SNP_pos), ]
BarEv <- POTevents.fun(T = -log10(org1$Pvalue), thres = 7, date = org1$SNP_pos)

# tB <-org1$SNP_pos/51304566
tB <-org1$SNP_pos
png("emplambda_eqtl_chr22_DrQ.png",height = 800,width = 1000)
emplambdaB <- emplambda.fun(posE = BarEv$Px,inddat=BarEv$inddat, t = tB, lint = 500, tit = "EUR_eqtl_chr22")
# dev.off()
# rs1<-cbind(Pi=BarEv$Pi,Px=BarEv$Px,datePx=BarEv$datePx,Im=BarEv$Im,Ix=BarEv$Ix,L=BarEv$L,inddat=BarEv$inddat,T=BarEv$T,chr=i)
rs1<-cbind(Px=BarEv$Px,datePx=BarEv$datePx,chr=i)
rs3<-cbind(inddat=BarEv$inddat,T=BarEv$T,chr=i)

#------------------------------------
rs2<-cbind(emplambda=emplambdaB$emplambda, t=emplambdaB$t,chr=i)
# output_file1<-paste("../output/NHPoisson_POTevents_Q_test.txt", sep = "")
# output_file2<-paste("../output/NHPoisson_emplambda_Q_test",j,".txt", sep = "")
write.table(rs1,"../output/NHPoisson_POTevents_Q_test_Px.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(rs3,"../output/NHPoisson_POTevents_Q_test_inddat.txt",row.names = F, col.names = T,quote =F,sep="\t")
# write.table(rs4,"../output/NHPoisson_POTevents_Q_test.txt",row.names = F, col.names = T,quote =F,sep="\t")

#---------------------
write.table(rs2,"../output/NHPoisson_emplambda_Q_test.txt",row.names = F, col.names = T,quote =F,sep="\t")

#------------------------------------the snp unique 


indx <-which(org$ QTL_type == "eQTL")
org1 <- org[indx, ][1:10000, ]

org_p <- org1%>%dplyr::select(SNP_pos,Pvalue)%>%unique()
org_pg <- group_by(org_p, SNP_pos)
delay<-summarise(org_pg,min_p = min(Pvalue))%>%as.data.frame()
org2<-delay
org2 <- org2[order(org2$SNP_pos), ]
BarEv <- POTevents.fun(T = -log10(org2$min_p), thres = 7, date = org2$SNP_pos)

tB <-org2$SNP_pos
png("emplambda_eqtl_chr22_DrQ_snp_unique.png",height = 800,width = 1000)
emplambdaB <- emplambda.fun(posE = BarEv$Px,inddat=BarEv$inddat, t = tB, lint = 500, tit = "EUR_eqtl_chr22_snp_unique")
dev.off()
rs1<-cbind(Pi=BarEv$Pi,Px=BarEv$Px,datePx=BarEv$datePx,Im=BarEv$Im,Ix=BarEv$Ix,L=BarEv$L,inddat=BarEv$inddat,T=BarEv$T,chr=i)
rs2<-cbind(emplambda=emplambdaB$emplambda, t=emplambdaB$t,chr=i)
# output_file1<-paste("../output/NHPoisson_POTevents_Q_test_snp_occur1",j,".txt", sep = "")
# output_file2<-paste("../output/NHPoisson_emplambda_Q_test_snp_occur1",j,".txt", sep = "")
write.table(rs1,"../output/NHPoisson_POTevents_Q_test_snp_unique.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(rs2,"../output/NHPoisson_emplambda_Q_test_snp_unique.txt",row.names = F, col.names = T,quote =F,sep="\t")