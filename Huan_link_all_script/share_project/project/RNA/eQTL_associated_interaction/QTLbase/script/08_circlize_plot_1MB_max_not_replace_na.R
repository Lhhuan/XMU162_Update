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
  #----------------------plot
  histogram_pwd<-paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",QTL,sep="")
  setwd(histogram_pwd)
  print(histogram_pwd)
  pdf("histogram_max_Emplambda_not_replace_na.pdf")
  p1<-ggplot(org_QTL, aes(x =max_emplambda)) +geom_histogram(position="identity",colour ="#222831")  #,fill ="black"
  p1<-p1+xlab("Emplambda") + ylab("Count")+p_theme
  p1<-p1+ggtitle(QTL) 
  p1<-p1+theme(plot.title = element_text(hjust = 0.5)) #title center

  print(p1)
  dev.off()
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
    rs <- bind_rows(rs,normalize_data)
    #-----------------------------------------
    pdf("histogram_max_Emplambda_not_replace_na_less.pdf")
    p1<-ggplot(AAAA1, aes(x =max_emplambda)) +geom_histogram(position="identity",colour ="#222831")  #,fill ="black"
    p1<-p1+xlab("Emplambda") + ylab("Count")+p_theme
    p1<-p1+ggtitle(QTL) 
    p1<-p1+theme(plot.title = element_text(hjust = 0.5)) #title center

    print(p1)
    dev.off()}
  else{
    normalize_data<- org_QTL[,c('chr','region','xQTL')]
    max_value = max(org_QTL$max_emplambda)
    min_value = min(org_QTL$max_emplambda)
    for(i in (1:nrow(org_QTL))){
      normalize_data[i,4]<-(org_QTL[i,4]-min_value)/(max_value - min_value) 
    }
    rs <- bind_rows(rs,normalize_data)}}


#------------------------
colnames(rs)[4] <- "max_emplambda"
library(stringr)
tmp <- str_split_fixed(rs$region,pattern='_',n=2)

rs$start<-tmp[,1]
rs$end<-tmp[,2]
data<- spread(rs,xQTL,max_emplambda)
normalize_data1<- data[,c("chr","start","end",'QTL','eQTL','mQTL','sQTL','caQTL','hQTL','pQTL','edQTL','reQTL','miQTL','cerQTL','metaQTL','riboQTL')]
#-------------
set.seed(999)
library("OmicCircos")
data("UCSC.hg19.chr")
#data("TCGA.BC.gene.exp.2k.60")
dat <- UCSC.hg19.chr
#dat$chrom <- gsub("chr", "",dat$chrom)
indx = which(dat$chrom != "chrX" & dat$chrom != "chrY")
seg.d = dat[indx, ]
seg.name = paste("chr", seq(1, 22, by = 1), sep = "")
seg.c = segAnglePo(seg.d, seg = seg.name)
# colors <- rainbow(10, alpha=0.5);
colors <- rainbow(22, alpha = 0.8)

#---------plot circos heatmap
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")
pdf(file = "emplambda_circos_heatmap_mean_interval_1MB_max_min_sacle_outlier_10_na_not_replace.pdf", height = 8, width = 8, compress = TRUE)
par(mar= c(1,1,1,1))
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "",main = "Interval = 1MB")

# circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos_refine(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE,col ='#2171b5')
color2<-c('#c6dbef','#6baed6','#2171b5')
# circos(R = 20, cir = seg.c, type= "heatmap", W = 380, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE,col ="red")
circos_refine(R = 100, cir = seg.c, type= "heatmap", W = 280, mapping = normalize_data1,B = T, huan_col=color2,col.v = 4,col.bar=TRUE)
dev.off()

#--------------