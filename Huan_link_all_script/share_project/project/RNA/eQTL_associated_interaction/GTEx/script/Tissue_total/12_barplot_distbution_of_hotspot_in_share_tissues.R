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
#-------------------------------------------------------
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_total/")

org<-read.table("11_share_tissue_number_count.txt.gz",header = T,sep = "\t") %>% as.data.frame()

color = "#5f939a"
pdf("barplot_distbution_of_hotspot_in_share_tissues.pdf",height = 3.5,width = 3.8)
p1 <-ggplot(org,mapping=aes(x=number_of_share_tissue,y=number_of_hotspots))+geom_bar(stat='identity',colour=color,fill =color,width=0.5 )+
    p_theme+xlab("Number of tissues") +ylab("Number of hotspots")
print(p1)
dev.off()

