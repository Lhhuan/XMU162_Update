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
org_a<-org_need%>%group_by(xQTL,chr,region)%>%summarise(max(emplambda))%>%as.data.frame()
colnames(org_a)[4] <- "max_emplambda"
org_b<- spread(org_a,xQTL,max_emplambda)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

QTLs<-unique(org_a$xQTL)
rs <- data.frame()
for(QTL in QTLs){
  print(QTL)
  org_QTL <- filter(org_a,xQTL == QTL)
  #------------------------------------------
  #max min normalize
  #------------------------------replace outlier 
  if(QTL=="QTL"|| QTL=="eQTL"||QTL=="mQTL"||QTL=="sQTL"||QTL=="caQTL"||QTL=="hQTL"||QTL=="pQTL"||QTL=="edQTL"){
    AAAA1<-org_QTL[order(org_QTL[,4],decreasing=T),] #SORT decreasing by QTL
    AAAA1[c(1:10),4] <-AAAA1[11,4] #----------------top 10 replace by 11
    normalize_data<- AAAA1[,c('chr','region','xQTL')]
    max_value = max(AAAA1$max_emplambda)
    min_value = min(AAAA1$max_emplambda)
    for(i in (1:nrow(AAAA1))){
      normalize_data[i,4]<-(AAAA1[i,4]-min_value)/(max_value - min_value) }
    rs <- bind_rows(rs,normalize_data)}
    #-----------------------------------------
  else{
    normalize_data<- org_QTL[,c('chr','region','xQTL')]
    max_value = max(org_QTL$max_emplambda)
    min_value = min(org_QTL$max_emplambda)
    for(i in (1:nrow(org_QTL))){
      normalize_data[i,4]<-(org_QTL[i,4]-min_value)/(max_value - min_value) 
    }
    rs <- bind_rows(rs,normalize_data)}}
colnames(rs)[4] <- "max_emplambda"
library(stringr)
tmp <- str_split_fixed(rs$region,pattern='_',n=2)

rs$start<-tmp[,1]
rs$end<-tmp[,2]
data<- spread(rs,xQTL,max_emplambda)
normalize_data1<- data[,c("chr","start","end",'QTL','eQTL','mQTL','sQTL','caQTL','hQTL','pQTL','edQTL','reQTL','miQTL','cerQTL','metaQTL','riboQTL')]

#----------------------------------------------------max min scale outlier  out
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/heatmap_chr/max_min_scale/less/")
for (i in c(1:22)){
  data_chr <-filter(data,chr ==i)
  data_chr_1<- data_chr[,c('region',"start","end",'QTL','eQTL','mQTL','sQTL','caQTL','hQTL','pQTL','edQTL','reQTL','miQTL','cerQTL','metaQTL','riboQTL')]
#--------------------------------------------------------------------------------------------------------------------------


  data_chr_1$start<-as.numeric(data_chr_1$start)
  data_chr_qqq<-data_chr_1[order(data_chr_1$start),]
  region<-data_chr_qqq$region

  data_chr_2<-data_chr_qqq[,-c(1,2,3)]
  colnames(data_chr_2)[1] <-"ALL QTL"
  row.names(data_chr_2)<-region

  data_used_heatmap<-as.matrix(data_chr_2)
  final_data<-t(data_used_heatmap)
  # color2<-c('#c6dbef','#6baed6','#2171b5')
  color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
  title= paste("Chr",i,sep =" ")
  figure_name<-paste("max_emplambda",i,"max_min_scale_outlier_out_not_replace_na.pdf",sep ="_")
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


