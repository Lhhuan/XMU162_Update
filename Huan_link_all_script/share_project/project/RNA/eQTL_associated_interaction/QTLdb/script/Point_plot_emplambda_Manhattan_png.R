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
    figure_name <-paste("NHPoisson_emplambda_interval_1000_cutoff_7.3_chr",i,".png", sep = "")
    title_name<-paste("Chr",i,sep = " ")
    png(figure_name,height = 1000,width = 2000)
    p1<-ggplot(org2, aes(x =t, y=emplambda,colour =xQTL)) +geom_point(size=0.1) 

    p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) 
    p2<-p2+xlab("Position") + ylab("Emplambda") 
    p2<-p2+ggtitle(title_name)
    p2<-p2+theme(legend.text =element_text(size =20),legend.title =element_text(size =20),plot.title=element_text(size =20))
    print(p2)
    dev.off()
    print(i)
}
