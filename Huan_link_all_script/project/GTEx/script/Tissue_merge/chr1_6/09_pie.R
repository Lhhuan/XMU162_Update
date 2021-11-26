# library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(parallel)
library(conflicted)
library(gridExtra)
library(Hmisc)

org<-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/communities.bed.gz",header =F,sep="\t")%>%data.frame()
colnames(org)[4] <-"Group"
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
a<-table(org$Group)%>%as.data.frame()


setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/")

pdf("6kmer_cluster_default.pdf")
pie(a$Freq, cex=1.5,col = mycolor,labels = a$Freq, radius = 1,main="Composition of hotspots",cex.main=2)
legend("topright", c("0","1","2","3","4","5"), cex = 1, fill = mycolor)
dev.off()
#---------------------------6
org<-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/communities_6.bed.gz",header =F,sep="\t")%>%data.frame()
colnames(org)[4] <-"Group"
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93","#1C7947")
a<-table(org$Group)%>%as.data.frame()


setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/")

pdf("6kmer_cluster6.pdf")
pie(a$Freq, cex=1.5,col = mycolor,labels = a$Freq, radius = 1,main="Composition of hotspots",cex.main=2)
legend("topright", c("0","1","2","3","4","5","6"), cex = 1, fill = mycolor)
dev.off()

#----------------------------5
org<-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/communities_5.bed.gz",header =F,sep="\t")%>%data.frame()
colnames(org)[4] <-"Group"
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
a<-table(org$Group)%>%as.data.frame()


setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/")

pdf("6kmer_cluster5.pdf")
pie(a$Freq, cex=1.5,col = mycolor,labels = a$Freq, radius = 1,main="Composition of hotspots",cex.main=2)
legend("topright", c("0","1","2","3","4","5"), cex = 1, fill = mycolor)
dev.off()

#----------------------------4
org<-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/communities_4.bed.gz",header =F,sep="\t")%>%data.frame()
colnames(org)[4] <-"Group"
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9")
a<-table(org$Group)%>%as.data.frame()


setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/")

pdf("6kmer_cluster4.pdf")
pie(a$Freq, cex=1.5,col = mycolor,labels = a$Freq, radius = 1,main="Composition of hotspots",cex.main=2)
legend("topright", c("0","1","2","3","4"), cex = 1, fill = mycolor)
dev.off()

#----------------------------4
org<-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/communities_3.bed.gz",header =F,sep="\t")%>%data.frame()
colnames(org)[4] <-"Group"
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF")
a<-table(org$Group)%>%as.data.frame()


setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/")

pdf("6kmer_cluster3.pdf")
pie(a$Freq, cex=1.5,col = mycolor,labels = a$Freq, radius = 1,main="Composition of hotspots",cex.main=2)
legend("topright", c("0","1","2","3"), cex = 1, fill = mycolor)
dev.off()

#-------------------------------------------------default5, sets 0

org<-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/communities_s0.bed.gz",header =F,sep="\t")%>%data.frame()
colnames(org)[4] <-"Group"
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
a<-table(org$Group)%>%as.data.frame()


setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/chr1_6/6kmer/")

pdf("6kmer_cluster_default_s0.pdf")
pie(a$Freq, cex=1.5,col = mycolor,labels = a$Freq, radius = 1,main="Composition of hotspots",cex.main=2)
legend("topright", c("0","1","2","3","4","5"), cex = 1, fill = mycolor)
dev.off()