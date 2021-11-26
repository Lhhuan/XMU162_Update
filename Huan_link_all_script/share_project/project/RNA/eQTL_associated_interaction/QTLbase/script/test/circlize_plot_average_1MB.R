load("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/OmicCircos_circos_refine")

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
org_a<-org_need%>%group_by(xQTL,chr,region)%>%summarise(mean(emplambda))%>%as.data.frame()
colnames(org_a)[4] <- "mean_emplambda"
data<- spread(org_a,xQTL,mean_emplambda)
nrow(data)
#3398

library(stringr)
tmp <- str_split_fixed(data$region,pattern='_',n=2)

data$start<-tmp[,1]
data$end<-tmp[,2]

data1<- data[,c("chr","start","end",'QTL','eQTL','mQTL','sQTL','caQTL','hQTL')]
data1[is.na(data1)]=0
data1$start<-as.integer(data1$start)
data1$end<-as.integer(data1$end)
# save(data1,file = "mean_emplambda_all")
#-------------------------------------------------
#plot circos using OmicCircos
library("OmicCircos")
#---------------------------------------------



#------------------------------------------------------ using UCSC hg19作为 background
# options(stringsAsFactors = FALSE)
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
#----------------------------------
n<-which.max(data1$caQTL)
data1<-data1[-n,]
scale_data1<- data1[,c('chr','start','end')]
normalize_data<- data1[,c('chr','start','end')]
#---------------------------scale data
for(i in c(4:ncol(data1))){
  for(j in c(1:nrow(data1))){
    max_value = max(data1[,i])
    min_value = min(data1[,i])
    times = 1/max_value
    scale_data1[,i]=times * data1[,i]
    # print(i)
    # print(times)
    normalize_data[j,i]<-(data1[j,i]-min_value)/(max_value - min_value)    }}

#--------------------------plot circos heatmap


pdf(file = "emplambda_circos_heatmap_mean_interval_1MB_times_scale.pdf", height = 8, width = 8, compress = TRUE)
# pdf(file = "test_circos_heatmap_new.pdf", height = 8, width = 8, compress = TRUE)
par(mar = c(1,1,1,1))
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "",main = "Interval = 1MB")

# circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos_refine(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE,col ='#2171b5' )
color2<-c('#c6dbef','#6baed6','#2171b5')
# circos(R = 20, cir = seg.c, type= "heatmap", W = 380, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE,col ="red")
circos_refine(R = 200, cir = seg.c, type= "heatmap", W = 180, mapping = normalize_data,B = T, huan_col=color2,col.v = 4,col.bar=TRUE)
dev.off()

pdf(file = "emplambda_circos_heatmap_mean_interval_1MB_max_min_scale.pdf", height = 8, width = 8, compress = TRUE)
# pdf(file = "test_circos_heatmap_new.pdf", height = 8, width = 8, compress = TRUE)
par(mar = c(1,1,1,1))
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "",main = "Interval = 1MB")

# circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos_refine(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE,col ='#2171b5' )
color2<-c('#c6dbef','#6baed6','#2171b5')
# circos(R = 20, cir = seg.c, type= "heatmap", W = 380, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE,col ="red")
circos_refine(R = 200, cir = seg.c, type= "heatmap", W = 180, mapping = normalize_data,B = T, huan_col=color2,col.v = 4,col.bar=TRUE)
dev.off()



#---------------------------------------------------

pdf(file = "emplambda_circos_heatmap_mean_interval_1MB.pdf", height = 8, width = 8, compress = TRUE)
par(mar = c(2,2,2,2))
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "",main = "Interval = 1MB")

# circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE)

circos(R = 200, cir = seg.c, type= "heatmap", W = 180, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE)
dev.off()
#-------------------------------------------------------------------------
# seg.f<-data1[,c("chr","start","end")]
# seg.f$the.v<-"NA"
# seg.f$NO<-"NA"
# seg.name =  seq(1, 22, by = 1)
# seg.c = segAnglePo(seg.f, seg = seg.name)
# colors <- rainbow(10, alpha=0.5);
# #-------------------------------------

# pdf(file = "emplambda_circos_heatmap_mean_interval_1MB.pdf", height = 8, width = 8, compress = TRUE)
# par(mar = c(2,2,2,2))
# plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "",main = "Interval = 1MB")

# # circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
# circos(R = 400, cir = seg.c, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE)

# circos(R = 20, cir = seg.c, type= "heatmap", W = 399, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE)
# dev.off()