library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)



setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/test_add_position/figure/")
org1<-fread("../output/NHPoisson_emplambda_interval_1000cutoff_7.3region.txt",header = T,sep = "\t") %>% as.data.frame()
org2<-fread("../output/NHPoisson_emplambda_interval_1000cutoff_7.3region_add.txt",header = T,sep = "\t") %>% as.data.frame()
org1$type<-c("QTLbase")
org2$type<-c("QTLbase_1KG")
final_file<-rbind(org1,org2)


png("emplambda_difference_between_org_and_add.png",height = 1000,width = 2000)
p1<-ggplot(final_file, aes(x =t, y=emplambda,colour =type)) +geom_point(size=0.1) 

p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),axis.line = element_line(colour = "black"),
    axis.text.y = element_text(size=20),axis.text.x = element_text(size=20)) 
p2<-p2+xlab("Position") + ylab("Emplambda") 
p2<-p2+ggtitle("Chr 22")
p2<-p2+theme(legend.text =element_text(size =20),legend.title =element_text(size =20),plot.title=element_text(size =20))
print(p2)
dev.off()
