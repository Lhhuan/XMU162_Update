library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(CMplot)



setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")

org<-read.table("../output/merge_QTL_all_QTLtype_pop.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# indx <-which(org$QTL_type == "eQTL")
EUR_eQTL <- filter(org,QTL_type == "eQTL" & Population == "EUR")
SNP <-paste(EUR_eQTL$SNP_chr,"_",EUR_eQTL$SNP_pos)
EUR_eQTL_data<-cbind(SNP,EUR_eQTL$SNP_chr,EUR_eQTL$SNP_pos,EUR_eQTL$Pvalue)

CMplot(EUR_eQTL_data,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=0.5,cir.legend=TRUE,outward=FALSE,cir.legend.col="black",cir.chr.h=0.6,chr.den.col="black",file="jpg",memo="EUR_eQTL_Manhattan_plot",dpi=300,file.output=TRUE,verbose=TRUE,threshold=7,LOG10=TRUE)

threshold=7
LOG10=TRUE
#-----------------------------
data(pig60K)
CMplot(pig60K,plot.type="c",chr.labels=paste("Chr",c(1:18,"X"),sep=""),r=0.4,cir.legend=TRUE,outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

