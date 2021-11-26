library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/chr1/factor_anno")
All <-read.table("../communities.bed.gz",header = F,sep = "\t") %>% as.data.frame()
CHROMATIN_Accessibility <-read.table("CHROMATIN_Accessibility.bed.gz",header = F,sep = "\t") %>% as.data.frame()
CTCF<-read.table("CTCF.bed.gz",header = F,sep = "\t") %>% as.data.frame()
TFBS<-read.table("TFBS.bed.gz",header = F,sep = "\t") %>% as.data.frame()

H3K27ac<-read.table("H3K27ac.bed.gz",header = F,sep = "\t") %>% as.data.frame()
H3K27me3<-read.table("H3K27me3.bed.gz",header = F,sep = "\t") %>% as.data.frame()
H3K36me3<-read.table("H3K36me3.bed.gz",header = F,sep = "\t") %>% as.data.frame()
H3K4me1<-read.table("H3K4me1.bed.gz",header = F,sep = "\t") %>% as.data.frame()
H3K4me3<-read.table("H3K4me3.bed.gz",header = F,sep = "\t") %>% as.data.frame()
H3K9ac<-read.table("H3K9ac.bed.gz",header = F,sep = "\t") %>% as.data.frame()
H3K9me3<-read.table("H3K9me3.bed.gz",header = F,sep = "\t") %>% as.data.frame()

# CHROMATIN_Accessibility$Clsss <-"CHROMATIN_Accessibility"
# CTCF$Clsss <-"CTCF"
# TFBS$Clsss <-"TFBS"
# H3K27ac$Clsss <-"H3K27ac"
# H3K27me3$Clsss <-"H3K27me3"
# H3K36me3$Clsss <-"H3K36me3"
# H3K4me1$Clsss <-"H3K4me1"
# H3K4me3$Clsss <-"H3K4me3"
# H3K9ac$Clsss <-"H3K9ac"
# H3K9me3$Clsss <-"H3K9me3"

#-----------------
i=0


a <-filter(All,V4==i)%>%dplyr::select(V1,V2,V3)
All_filter_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(CTCF,V4==i)%>%dplyr::select(V1,V2,V3)
CTCF_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(TFBS,V4==i)%>%dplyr::select(V1,V2,V3)
TFBS_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(H3K27ac,V4==i)%>%dplyr::select(V1,V2,V3)
H3K27ac_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(H3K27me3,V4==i)%>%dplyr::select(V1,V2,V3)
H3K27me3_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(H3K36me3,V4==i)%>%dplyr::select(V1,V2,V3)
H3K36me3_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(H3K4me1,V4==i)%>%dplyr::select(V1,V2,V3)
H3K4me1_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(H3K4me3,V4==i)%>%dplyr::select(V1,V2,V3)
H3K4me3_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(H3K9ac,V4==i)%>%dplyr::select(V1,V2,V3)
H3K9ac_key <-paste(a$V1,a$V2,a$V3,sep="_")

a <-filter(H3K9me3,V4==i)%>%dplyr::select(V1,V2,V3)
H3K9me3_key <-paste(a$V1,a$V2,a$V3,sep="_")


setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/chr1/figure/")

library(VennDiagram)

T<-venn.diagram(list(All=All_filter_key,CTCF=CTCF_key,TFBS=TFBS_key),filename=NULL,lwd=1,lty=2,col=c('red','green','blue') ,fill=c('red','green','blue'),cat.col=c('red','green','blue'),reverse=TRUE)

grid.draw(T)



T<-venn.diagram(list(All=All_filter_key,CTCF=CTCF_key,TFBS=TFBS_key,H3K27ac=H3K27ac_key,H3K27me3=H3K27me3_key,H3K36me3=H3K36me3_key,H3K4me1=H3K4me1_key,H3K4me3=H3K4me3_key),filename=NULL,lwd=1,lty=2,col=c('yellow','red','red','red','red','red','green','blue') ,fill=c('yellow','red','red','red','red','red','green','blue'),cat.col=c('yellow','red','red','red','red','red','green','blue'),reverse=TRUE)

grid.draw(T)
#-----------------
T<-venn.diagram(list(All=All_filter_key,CTCF=CTCF_key,TFBS=TFBS_key,H3K27ac=H3K27ac_key,H3K27me3=H3K27me3_key,H3K36me3=H3K36me3_key,H3K4me1=H3K4me1_key,H3K4me3=H3K4me3_key),filename=NULL,lwd=1,lty=2,col=c('yellow','red','caf7e3','edffec','f6dfeb','e4bad4','green','blue') ,fill=c('yellow','red','caf7e3','edffec','f6dfeb','e4bad4','green','blue'),cat.col=c('yellow','red','caf7e3','edffec','f6dfeb','e4bad4','green','blue'),reverse=TRUE)

grid.draw(T)

T<-venn.diagram(list(All=All_filter_key,CTCF=CTCF_key,TFBS=TFBS_key,H3K27ac=H3K27ac_key,H3K27me3=H3K27me3_key),filename=NULL,lwd=1,lty=2,col=c('yellow','red','#bb8082','blue','#b8b5ff') ,fill=c('yellow','red','#bb8082','blue','#b8b5ff'),cat.col=c('yellow','red','#bb8082','blue','#b8b5ff'),reverse=TRUE)

grid.draw(T)
dev.off()


T<-venn.diagram(list(CTCF=CTCF_key,TFBS=TFBS_key,H3K27ac=H3K27ac_key,H3K27me3=H3K27me3_key),filename=NULL,lwd=1,lty=2,col=c('red','#bb8082','blue','#b8b5ff') ,fill=c('red','#bb8082','blue','#b8b5ff'),cat.col=c('red','#bb8082','blue','#b8b5ff'),reverse=TRUE)

grid.draw(T)
dev.off()