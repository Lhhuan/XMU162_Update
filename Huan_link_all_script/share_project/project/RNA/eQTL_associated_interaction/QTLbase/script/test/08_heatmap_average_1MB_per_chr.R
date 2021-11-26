source("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/OmicCircos_circos_refine_revise_color.R")

#-------------------------------refine function to revise the color of heatmap and heatmap2

library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(circlize)
library(CMplot)
library(reshape2)
library(tidyverse)
library(pheatmap)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")
org<-fread("../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_1MB_bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(org)[1] <- "chr"
colnames(org)[2] <- "start1"
colnames(org)[3] <- "end1"
colnames(org)[5] <- "start2"
colnames(org)[6] <- "end2"
colnames(org)[7] <- "xQTL"
colnames(org)[8] <- "emplambda"
org_need<-org[,c('chr','start1','end1','start2','end2','xQTL','emplambda')]
org_need$region<-paste(org_need$start1,org_need$end1,sep="_")
org_a<-org_need%>%group_by(xQTL,chr,region)%>%summarise(mean(emplambda))%>%as.data.frame()
colnames(org_a)[4] <- "mean_emplambda"


data<- spread(org_a,xQTL,mean_emplambda)
nrow(data)
#2733
data[is.na(data)]=0
colnames(data)[3] <-"ALL QTL"
#---------------------------------------------original
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/heatmap_chr/orginal/")
for (i in c(1:22)){
  data_chr <-filter(data,chr ==i)
  data_chr_1<- data_chr[,c("region",'ALL QTL','eQTL','mQTL','sQTL','caQTL','hQTL')]

  tmp <- str_split_fixed(data_chr_1$region,pattern='_',n=2)

  data_chr_1$start<-tmp[,1]
  data_chr_1$end<-tmp[,2]
  data_chr_1$start<-as.numeric(data_chr_1$start)
  data_chr_qqq<-data_chr_1[order(data_chr_1$start),]
  region<-data_chr_qqq$region

  data_chr_2<-data_chr_qqq[,-c(1,8,9)]
  row.names(data_chr_2)<-region

  data_used_heatmap<-as.matrix(data_chr_2)
  final_data<-t(data_used_heatmap)
  # color2<-c('#c6dbef','#6baed6','#2171b5')
  color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
  title= paste("Chr",i,sep =" ")
  figure_name<-paste("Chr",i,"orginal.pdf",sep ="_")
  pdf(figure_name)
  pheatmap(data_used_heatmap,cluster_rows = FALSE, cluster_cols = FALSE,color = color2,fontsize_row=1.5,fontsize_col=4,cellwidth=7,fontsize=5, main =title)
  dev.off()}

#-----------------------------------------------------times_scale and max min scale original
data1<- data[,c("chr","region",'ALL QTL','eQTL','mQTL','sQTL','caQTL','hQTL')]



scale_data1<- data1[,c("chr","region")]
normalize_data<- data1[,c("chr","region")]

for(i in c(3:ncol(data1))){
  for(j in c(1:nrow(data1))){
    max_value = max(data1[,i])
    min_value = min(data1[,i])
    times = 1/max_value
    scale_data1[,i]=times * data1[,i]
    # print(i)
    # print(times)
    normalize_data[j,i]<-(data1[j,i]-min_value)/(max_value - min_value)    }}

colnames(scale_data1)[3] <- 'ALL QTL'
colnames(scale_data1)[4] <- 'eQTL'
colnames(scale_data1)[5] <- 'mQTL'
colnames(scale_data1)[6] <- 'sQTL'
colnames(scale_data1)[7] <- 'caQTL'
colnames(scale_data1)[8] <- 'hQTL'
#---------------------------------
colnames(normalize_data)[3] <- 'ALL QTL'
colnames(normalize_data)[4] <- 'eQTL'
colnames(normalize_data)[5] <- 'mQTL'
colnames(normalize_data)[6] <- 'sQTL'
colnames(normalize_data)[7] <- 'caQTL'
colnames(normalize_data)[8] <- 'hQTL'
#-------------------times_scale original
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/heatmap_chr/times_scale/original/")
for (i in c(1:22)){
  data_chr <-filter(scale_data1,chr ==i)
  data_chr_1<- data_chr[,c("region",'ALL QTL','eQTL','mQTL','sQTL','caQTL','hQTL')]

  tmp <- str_split_fixed(data_chr_1$region,pattern='_',n=2)

  data_chr_1$start<-tmp[,1]
  data_chr_1$end<-tmp[,2]
  data_chr_1$start<-as.numeric(data_chr_1$start)
  data_chr_qqq<-data_chr_1[order(data_chr_1$start),]
  region<-data_chr_qqq$region

  data_chr_2<-data_chr_qqq[,-c(1,8,9)]
  row.names(data_chr_2)<-region

  data_used_heatmap<-as.matrix(data_chr_2)
  final_data<-t(data_used_heatmap)
  # color2<-c('#c6dbef','#6baed6','#2171b5')
  color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
  title= paste("Chr",i,sep =" ")
  figure_name<-paste("Chr",i,"times_scale_original.pdf",sep ="_")
  pdf(figure_name)
  pheatmap(data_used_heatmap,cluster_rows = FALSE, cluster_cols = FALSE,color = color2,fontsize_row=1.5,fontsize_col=4,cellwidth=7,fontsize=5, main =title)
  dev.off()}
#----------------------------------------------------max_min_scale original
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/heatmap_chr/max_min_scale/original/")
for (i in c(1:22)){
  data_chr <-filter(normalize_data,chr ==i)
  data_chr_1<- data_chr[,c("region",'ALL QTL','eQTL','mQTL','sQTL','caQTL','hQTL')]

  tmp <- str_split_fixed(data_chr_1$region,pattern='_',n=2)

  data_chr_1$start<-tmp[,1]
  data_chr_1$end<-tmp[,2]
  data_chr_1$start<-as.numeric(data_chr_1$start)
  data_chr_qqq<-data_chr_1[order(data_chr_1$start),]
  region<-data_chr_qqq$region

  data_chr_2<-data_chr_qqq[,-c(1,8,9)]
  row.names(data_chr_2)<-region

  data_used_heatmap<-as.matrix(data_chr_2)
  final_data<-t(data_used_heatmap)
  # color2<-c('#c6dbef','#6baed6','#2171b5')
  color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
  title= paste("Chr",i,sep =" ")
  figure_name<-paste("Chr",i,"max_min_scale_original.pdf",sep ="_")
  pdf(figure_name)
  pheatmap(data_used_heatmap,cluster_rows = FALSE, cluster_cols = FALSE,color = color2,fontsize_row=1.5,fontsize_col=4,cellwidth=7,fontsize=5, main =title)
  dev.off()}

#----------------------------------------------------times_scale and max min scale less
AAAA<-data1
colnames(AAAA)[3] <- 'QTL'
AAAA1<-AAAA[order(AAAA$QTL,decreasing=T),] #SORT decreasing by QTL
AAAA1 <-AAAA1[-c(1:4),]    #remove the top 3 row
AAAA2 <-AAAA1[order(AAAA1$sQTL,decreasing=T),]  
AAAA2 <-AAAA2[-c(1:4),]
AAAA3 <-AAAA2[order(AAAA2$caQTL,decreasing=T),]
AAAA3 <-AAAA3[-c(1:13),]
AAAA4 <-AAAA3[order(AAAA3$hQTL,decreasing=T),]
AAAA4 <-AAAA4[-c(1:4),]
data1<-AAAA4
#-------------------------------
scale_data1<- data1[,c("chr","region")]
normalize_data<- data1[,c("chr","region")]

for(i in c(3:ncol(data1))){
  for(j in c(1:nrow(data1))){
    max_value = max(data1[,i])
    min_value = min(data1[,i])
    times = 1/max_value
    scale_data1[,i]=times * data1[,i]
    # print(i)
    # print(times)
    normalize_data[j,i]<-(data1[j,i]-min_value)/(max_value - min_value)    }}

colnames(scale_data1)[3] <- 'ALL QTL'
colnames(scale_data1)[4] <- 'eQTL'
colnames(scale_data1)[5] <- 'mQTL'
colnames(scale_data1)[6] <- 'sQTL'
colnames(scale_data1)[7] <- 'caQTL'
colnames(scale_data1)[8] <- 'hQTL'
#---------------------------------
colnames(normalize_data)[3] <- 'ALL QTL'
colnames(normalize_data)[4] <- 'eQTL'
colnames(normalize_data)[5] <- 'mQTL'
colnames(normalize_data)[6] <- 'sQTL'
colnames(normalize_data)[7] <- 'caQTL'
colnames(normalize_data)[8] <- 'hQTL'
#-------------------times_scale less
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/heatmap_chr/times_scale/less/")
for (i in c(1:22)){
  data_chr <-filter(scale_data1,chr ==i)
  data_chr_1<- data_chr[,c("region",'ALL QTL','eQTL','mQTL','sQTL','caQTL','hQTL')]

  tmp <- str_split_fixed(data_chr_1$region,pattern='_',n=2)

  data_chr_1$start<-tmp[,1]
  data_chr_1$end<-tmp[,2]
  data_chr_1$start<-as.numeric(data_chr_1$start)
  data_chr_qqq<-data_chr_1[order(data_chr_1$start),]
  region<-data_chr_qqq$region

  data_chr_2<-data_chr_qqq[,-c(1,8,9)]
  row.names(data_chr_2)<-region

  data_used_heatmap<-as.matrix(data_chr_2)
  final_data<-t(data_used_heatmap)
  # color2<-c('#c6dbef','#6baed6','#2171b5')
  color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
  title= paste("Chr",i,sep =" ")
  figure_name<-paste("Chr",i,"times_scale_less.pdf",sep ="_")
  pdf(figure_name)
  pheatmap(data_used_heatmap,cluster_rows = FALSE, cluster_cols = FALSE,color = color2,fontsize_row=1.5,fontsize_col=4,cellwidth=7,fontsize=5, main =title)
  dev.off()}
#----------------------------------------------------max min scale less
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/heatmap_chr/max_min_scale/less/")
for (i in c(1:22)){
  data_chr <-filter(normalize_data,chr ==i)
  data_chr_1<- data_chr[,c("region",'ALL QTL','eQTL','mQTL','sQTL','caQTL','hQTL')]

  tmp <- str_split_fixed(data_chr_1$region,pattern='_',n=2)

  data_chr_1$start<-tmp[,1]
  data_chr_1$end<-tmp[,2]
  data_chr_1$start<-as.numeric(data_chr_1$start)
  data_chr_qqq<-data_chr_1[order(data_chr_1$start),]
  region<-data_chr_qqq$region

  data_chr_2<-data_chr_qqq[,-c(1,8,9)]
  row.names(data_chr_2)<-region

  data_used_heatmap<-as.matrix(data_chr_2)
  final_data<-t(data_used_heatmap)
  # color2<-c('#c6dbef','#6baed6','#2171b5')
  color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
  title= paste("Chr",i,sep =" ")
  figure_name<-paste("Chr",i,"max_min_scale_less.pdf",sep ="_")
  pdf(figure_name)
  pheatmap(data_used_heatmap,cluster_rows = FALSE, cluster_cols = FALSE,color = color2,fontsize_row=1.5,fontsize_col=4,cellwidth=7,fontsize=5, main =title)
  dev.off()}

#-------------------------------------------




#---------------------------------------------
data2<-data1[order(data1$start),]

tmp <- str_split_fixed(data1$region,pattern='_',n=2)

data1$start<-tmp[,1]
data1$end<-tmp[,2]
data1$start<-as.numeric(data1$start)










#---------------------------------------ggplot for heatmap

aa<-filter(org_a ,chr==22)
p <- ggplot(aa, aes(x=region,y=xQTL))+ geom_tile(aes(fill=mean_emplambda)) 
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1, vjust=1))
p <- p + scale_fill_gradient(low = '#c6dbef', high = '#2171b5')
p<-p+scale_y_discrete(breaks=c('ALL_QTL','eQTL','mQTL','sQTL','caQTL','hQTL','cerQTL','edQTL','riboQTL','reQTL','pQTL','miQTL','metaQTL','lncRNAQTL')
print(p)
dev.off()


