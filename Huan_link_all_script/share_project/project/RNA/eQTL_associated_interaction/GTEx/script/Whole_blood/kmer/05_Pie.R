library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)
library(gridExtra)
library(Hmisc)

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/")
org<-read.table("communities.bed.gz",header =F,sep="\t")%>%data.frame()
colnames(org)[4] <-"Group"
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
a<-table(org$Group)%>%as.data.frame()


pdf("./figure/kmer_cluster.pdf")
pie(a$Freq, cex=1.5,col = mycolor,labels = a$Freq, radius = 1,main="Composition of hotspots",cex.main=2)
legend("topright", c("0","1","2","3","4","5"), cex = 1, fill = mycolor)
dev.off()
