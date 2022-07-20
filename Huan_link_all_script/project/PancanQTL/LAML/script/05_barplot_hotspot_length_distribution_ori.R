[library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(tidyverse)
library(reshape2)
library(Seurat)
library(circlize)


p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
                                                
setwd("/home/huanhuan/project/PancanQTL/LAML/output/figure/")
org<-read.table("../04_LAML_cis_eQTL_hotspot.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)<-c("chr","start","end")
org$hotspot_length <-org$end - org$start
org$Log10_len <-log10(org$hotspot_length)

#---------------------------

# a <-as.data.frame(table(org$Log10_len))
# a<-a[order(a$Freq,a$Var1),]
# a$Var1 <- as.numeric(as.character(a$Var1))
# a$Freq <- as.numeric(as.character(a$Freq))
# llen_c <-a[38,"Var1"]
#3.35983548233989
#----------------------------
org$class <-NA
org$class[which(org$hotspot_length ==1 )] <- "1"
org$class[which(org$hotspot_length <=10 & org$hotspot_length >1)] <- "2-10"
org$class[which(org$hotspot_length <=100 & org$hotspot_length >10)] <- "11-100"
org$class[which(org$hotspot_length <=500 & org$hotspot_length >100)] <- "101-500"
org$class[which(org$hotspot_length <=1000 & org$hotspot_length >500)] <- "501-1,000"
org$class[which(org$hotspot_length <=1500 & org$hotspot_length >1000)] <- "1,001-1,500"
org$class[which(org$hotspot_length <=2000 & org$hotspot_length >1500)] <- "1,501-1670"


org$class<-factor(org$class,levels=c("1","2-10","11-100","101-500","501-1,000","1,001-1,500","1,501-1670"))
l_c <- table(org$class)%>%data.frame()
colnames(l_c) <-c("Length","Frequency")

write.table(l_c,"../05_ori_hotspot_distribution_bar.txt",row.names=F,quote=F,sep="\t")

pdf("./05_ori_barplot_hotspot_distribution.pdf",width = 4,height = 4)
p1 <-ggplot(data = l_c, mapping = aes(x =Length , y = Frequency)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.6)+
  p_theme+scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black",angle=60,hjust=1))
print(p1)
dev.off()
#------------------------------filter top length
# filter_org <-filter(org,Log10_len < as.numeric(as.character(llen_c)))


# write.table(filter_org[,1:3],"../../../output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed",row.names=F,col.names =F,quote=F,sep="\t")


]