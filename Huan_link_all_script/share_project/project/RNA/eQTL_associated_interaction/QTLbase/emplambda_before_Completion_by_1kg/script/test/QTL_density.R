library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")

org<-read.table("../output/EUR_Blood_mQTL.txt",header = T,sep = "\t") %>% as.data.frame()

p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))


pdf("Density_of_QTL_region.pdf",height = 4,width = 5) 
p1 <-ggplot (org,aes(x=position,colour = as.factor(chr)))+ #不把chr转成factor就不会分染色体，因为chr那列是纯数字
    geom_density(alpha=.2) + 
    xlab("Position")+
    p.theme+
    ylab("Density")+
    labs(colour = "Chr") ##修改图例名字
print(p1)
dev.off()
