library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tidyverse)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_total/")

all_hotspot <- read.table("number_of_tissue_all_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()
tissue_specific_hotspot <- read.table("number_of_tissue_specific_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()


all_data<-inner_join(all_hotspot,tissue_specific_hotspot,by="Tissue")
all_data$diff <- all_data$Number.x - all_data$Number.y
diff <-all_data%>%dplyr::select(-c(2,3))
colnames(diff)[2] <-"Number"

tissue_specific_hotspot$Class <- "Tissue specific hotspot"
diff$Class <- "All hotspot"
final <-bind_rows(tissue_specific_hotspot,diff)
final$Tissue<-factor(final$Tissue,levels= sort(diff$Tissue,decreasing = T))
pdf("13_tissue_specific_and_all_tissue.pdf")

# pdf("13_tissue_specific_and_all_tissue.pdf",height = 7,width = 7.5)
p1 <-ggplot(final,aes(x=Tissue,y=Number,fill=Class))+
    geom_bar(stat='identity',width=0.5 )+
    scale_fill_manual(values=c("#f29191","#1eae98"))+
    p_theme+xlab("Tissue") +ylab("Number of hotspots") +
    coord_flip()+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=7))

print(p1)   
dev.off()


all_hotspot$Tissue<-factor(all_hotspot$Tissue,levels= sort(all_hotspot$Tissue,decreasing = T))
pdf("13_number_of_hotspot_tissue.pdf")

# pdf("13_tissue_specific_and_all_tissue.pdf",height = 7,width = 7.5)
p1 <-ggplot(all_hotspot,aes(x=Tissue,y=Number,fill="#1eae98"))+
    geom_bar(stat='identity',width=0.5)+
    # scale_fill_manual(values="#1eae98")+
    p_theme+xlab("Tissue") +ylab("Number of hotspots") +
    coord_flip()+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"),
   legend.position = "none" )

print(p1)   
dev.off()
