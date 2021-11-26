library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))


xQTL = "eQTL"
cutoff= 7.3
tissue = "Whole_Blood"
# intervals=c(6,7,8,9,12,15,18)
j=18


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/NHP/")
org<-read.table(paste0("NHPoisson_emplambda_interval_",j,"_cutoff_",cutoff,"_",tissue,".txt.gz"),header = T,sep = "\t") %>% as.data.frame()

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/figure/Whole_Blood/")
pdf("04_NHPoisson_emplambda_violin.pdf",width=3.5, height=3.5)
p<-ggplot(org,aes(x=1, y=emplambda))+geom_violin(fill="#a3d2ca",width=0.65)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("emplambda")

print(p)
dev.off()


+geom_boxplot(fill = "#5eaaa8",width=0.2)