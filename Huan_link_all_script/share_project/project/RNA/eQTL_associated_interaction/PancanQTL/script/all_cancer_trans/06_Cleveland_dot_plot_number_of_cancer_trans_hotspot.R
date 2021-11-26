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
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/trans/")
org<-read.table("05_number_of_cancer_all_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org_filter <- org%>%filter(Number >10)

# color = "#5f939a"
color = "#008c9e"
p2 <-ggplot(org_filter,mapping=aes(x=cancer,y=Number)) +
    p_theme+geom_segment(aes(xend =cancer),yend = 0,stat = "identity",colour="#686868") +
    geom_point(colour=color)+
    xlab("")+
    theme(axis.text.y = element_text(size =6,color="black"),
    axis.title.y = element_text(size =6,color="black"),
    axis.text.x = element_text(size =5.5,color="black",angle = 30,hjust=1),
    legend.title=element_blank(),
    legend.text=element_text(size=8),
    legend.key.size=unit(0.3, "lines"))

pdf("06_Cleveland_dot_plot_number_of_cancer_trans_hotspot.pdf",height = 3,width = 4)
print(p2)
dev.off()

