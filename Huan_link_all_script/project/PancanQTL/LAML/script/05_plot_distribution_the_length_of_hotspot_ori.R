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
setwd("/home/huanhuan/project/PancanQTL/LAML/output/figure/")
org<-read.table("../04_LAML_cis_eQTL_hotspot.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)[1] <-"CHR"
colnames(org)[2] <-"start"
colnames(org)[3] <-"end"
org$hotspot_length <-org$end - org$start

pdf("05_ori_boxplot_distribution_extend_hotspot_length_log10.pdf",width=3.5, height=3.5)
p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("log10(length of hotspot)") #+coord_cartesian(ylim = c(0, 6.5))

print(p)
dev.off()

pdf("05_ori_boxplot_distribution_extend_hotspot_length_log10_no_ylim.pdf",width=3.5, height=3.5)
p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("log10(length of hotspot)")

print(p)
dev.off()

stat <-data.frame(mean = mean(org$hotspot_length),median=median(org$hotspot_length),min=min(org$hotspot_length),max=max(org$hotspot_length),sd=sd(org$hotspot_length))%>%t()%>%data.frame()
colnames(stat) <-"length"
stat <-add_column(stat, statistic=rownames(stat), .before = "length")

write.table(stat,"../05_ori_distribution_of_extend_hotspot_length.txt",row.names = F, col.names = T,quote =F,sep="\t")
