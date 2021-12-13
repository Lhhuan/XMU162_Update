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
org<-read.table("../../../output/Tissue_merge/Cis_eQTL/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833_GC.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)[1:5] <-c("chr","start","end","AT_content","GC_content")
pdf("031_violin_GC_content.pdf",width=3.5, height=3.5)
p<-ggplot(org,aes(x=1, y=GC_content))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("GC content") #+coord_cartesian(ylim = c(0, 6.5))

print(p)
dev.off()

