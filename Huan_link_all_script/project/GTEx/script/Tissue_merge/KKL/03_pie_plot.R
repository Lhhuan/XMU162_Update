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

org<-read.table("/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/01_kkl_result.txt",header =T,sep="\t")%>%data.frame()
# colnames(org)[4] <-"Group"
mycolor<-c( "#70A1D7","#C86B85","#FFD2A5", "#C79ECF", "#C1C0B9", "#A1DE93")
a<-table(org$cluster)%>%as.data.frame()


setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/KKL/")

pdf("kkl_cluster_default.pdf")
pie(a$Freq, cex=1.5,col = mycolor,labels = a$Freq, radius = 1,main="Composition of hotspots",cex.main=2)
legend("topright", c("1","2","3","4","5","6"), cex = 1, fill = mycolor)
dev.off()