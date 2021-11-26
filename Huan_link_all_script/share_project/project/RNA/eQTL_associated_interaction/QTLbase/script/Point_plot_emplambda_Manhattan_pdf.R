library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)



setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/xQTL_merge_Manhattan/")
org<-fread("../../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()

for(i in c(1:22)){
    org2 <-filter(org,chr==i)
    org2$emplambda<-as.numeric(org2$emplambda)
    org2$t<-as.numeric(org2$t)
    figure_name <-paste("NHPoisson_emplambda_interval_1000_cutoff_7.3_chr",i,".pdf", sep = "")
    title_name<-paste("Chr",i,sep = " ")
    pdf(figure_name,height = 5,width = 10)
    p1<-ggplot(org2, aes(x =t, y=emplambda,colour =xQTL)) +geom_point(size=0.001) 

    p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10),axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10),axis.text.x = element_text(size=10)) 
    p2<-p2+xlab("Position") + ylab("Emplambda") 
    p2<-p2+ggtitle(title_name)
    print(p2)
    dev.off()
    print(i)
}