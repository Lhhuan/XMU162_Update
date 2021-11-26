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
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")

org<-read.table("11_count_cancer_tissue_share_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$cancer_tissue <-paste(org$cancer,org$tissue,sep=":")

org$cancer_tissue<-factor(org$cancer_tissue,levels= sort(org$cancer_tissue,decreasing = T))

p1 <- ggplot(org,mapping=aes(x=cancer_tissue,y=Number)) +
    geom_bar(fill= "#008c9e",stat = "identity")+
    ylab("Number of hotspots")+
    coord_flip() +p_theme



pdf("12_barplot_number_of_cancer_share_hotspot.pdf")
print(p1)
dev.off()







p2 <-ggplot(org,mapping=aes(x=cancer_tissue,y=Number)) +
    p_theme+geom_segment(aes(xend =cancer_tissue),yend = 0,colour = 'grey35')+
    geom_point(size = 2, color= "#008c9e")+
    xlab("")+
    theme(axis.text.y = element_text(size =9,color="black"),
    axis.text.x = element_text(size =5.4,color="black",angle = 30,hjust=1))
    # scale_y_continuous(expand = c(0,0))
pdf("12_Cleveland_dot_plot_number_of_cancer_share_hotspot.pdf",height = 3,width = 5.5)
print(p2)
dev.off()
