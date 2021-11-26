library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)

setwd("~/project/RNA/eQTL_associated_interaction/QTLbase/figure/")
QTLs = c("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL");
for (xQTL in QTLs){
    file<-paste("../output/ALL_",xQTL,"/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",xQTL,"_cis_trans_distance.txt.gz",sep="")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    figure_name <- paste("./ALL_",xQTL,"/distance_emplambda_point.pdf",sep="")
    pdf(figure_name,height = 7,width = 7)
    p1<-ggplot(org, aes(x =distance, y=emplambda)) +geom_point(size=0.01, color ="black" ) 

    p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10),axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10),axis.text.x = element_text(size=10)) 
    p2<-p2+xlab("Distance") + ylab("Emplambda") 
    p2<-p2+ggtitle(xQTL)
    print(p2)
    dev.off()
    print(xQTL)

}


for (xQTL in QTLs){
    file<-paste("../output/ALL_",xQTL,"/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_",xQTL,"_cis_trans_distance.txt.gz",sep="")
    org<-fread(file,header = T,sep = "\t") %>% as.data.frame()
    figure_name <- paste("./ALL_",xQTL,"/distance_emplambda_point.png",sep="")
    png(figure_name,height = 1000,width = 1000)
    p1<-ggplot(org, aes(x =distance, y=emplambda)) +geom_point(size=0.1, color ="black" ) 

    p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) 
    p2<-p2+xlab("Distance") + ylab("Emplambda") 
    p2<-p2+ggtitle(xQTL)
    print(p2)
    dev.off()
    print(xQTL)
}
