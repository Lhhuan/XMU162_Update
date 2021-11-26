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

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

#---------------------------------------------------------
library(Hmisc)

org<-read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org) <-c("CHR","start","end")
org$hotspot_length <-org$end - org$start
org$Log10_len <-log10(org$hotspot_length)
org$chr=gsub("chr","",org$CHR)
# i=1
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")
plot_log10 <-function(i){
    rs <-filter(org,chr==i)
    p<-ggplot(rs,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ 
        theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(size=5,color="black"),plot.title = element_text(hjust = 0.5))+xlab("")+ylab("log10(length of hotspot)")+
        coord_cartesian(ylim = c(0, 6.5))+scale_y_continuous(breaks=seq(0,6.5,1))+ggtitle(paste0("Chr",i)) 
}

plist = lapply(1:22,plot_log10)
pdf("0411_log10_length_distribution_chr.pdf",width=10.3, height=8)
CombinePlots(plist,ncol=6,nrow=4)
dev.off()
