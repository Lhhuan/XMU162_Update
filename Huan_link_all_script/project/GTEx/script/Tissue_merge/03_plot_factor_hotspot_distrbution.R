library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(Seurat)
library(R.utils)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

#---------------------------------------------------------
library(Hmisc)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")
org<-read.table("../../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$chr,org$start,org$end,sep="_")
org <-org[,c("Marker","overlap_fraction","hotspot")]
org1 <-org%>%group_by(Marker,hotspot)%>%summarise(sum_overlap_fraction =sum(overlap_fraction))%>%data.frame()

org1$Marker <-gsub("CHROMATIN_Accessibility","CA",org1$Marker)
mk <- c("CA","TFBS","CTCF","H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3")



# pdf("03_violin_factor_fraction_of_hotspot_CA.pdf",width=5, height=5)
plot <-function(marker=NULL){
    dat <-filter(org1,Marker == marker)
    p<-ggplot(dat,aes(x=Marker, y=sum_overlap_fraction))+geom_violin(fill="#90A0CF",width=1)+geom_boxplot(fill = "#8BB6E0",width=0.04,outlier.color=NA)+ theme(legend.position ="none")+p_theme+xlab("")+ylab("Fraction")+ggtitle(marker)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5))
}
plist <-lapply(mk,plot)
pdf("./03_violin_boxplot_factor_fraction_of_hotspot.pdf",width=8.6, height=6)
CombinePlots(plist,ncol=4,nrow=3)
dev.off()

plot <-function(marker=NULL){
    dat <-filter(org1,Marker == marker)
    p<-ggplot(dat,aes(x=Marker, y=sum_overlap_fraction))+geom_boxplot(fill = "#8BB6E0",width= 0.5,outlier.color=NA)+ theme(legend.position ="none")+p_theme+xlab("")+ylab("Fraction")+ggtitle(marker)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5))
}
plist <-lapply(mk,plot)
pdf("./03_boxplot_factor_fraction_of_hotspot.pdf",width=8.6, height=6)
CombinePlots(plist,ncol=4,nrow=3)
dev.off()

plot <-function(marker=NULL){
    dat <-filter(org1,Marker == marker)
    p<-ggplot(dat,aes(x=Marker, y=sum_overlap_fraction))+geom_violin(fill = "#90A0CF",width=0.5)+ theme(legend.position ="none")+p_theme+xlab("")+ylab("Fraction")+ggtitle(marker)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),plot.title = element_text(hjust = 0.5))
}
plist <-lapply(mk,plot)
pdf("./03_violin_factor_fraction_of_hotspot.pdf",width=8.6, height=6)
CombinePlots(plist,ncol=4,nrow=3)
dev.off()