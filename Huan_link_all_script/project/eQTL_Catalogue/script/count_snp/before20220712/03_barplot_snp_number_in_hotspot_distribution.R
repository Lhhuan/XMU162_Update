library(ggplot2)
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


p_theme<-theme(panel.grid =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.line = element_line(colour = "black"))
                                 
setwd("/home/huanhuan/project/eQTL_Catalogue/script/count_snp/output/")
org<-read.table("02_count_snp_in_hotspot.bed.gz",header = T,sep = "\t") %>% as.data.frame()



#----------------------------
org$class <-org$SNP_number
org$class[which(org$SNP_number >=18 )] <- ">=18"


org$class<-factor(org$class,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17",">=18"))
l_c <- table(org$class)%>%data.frame()
colnames(l_c) <-c("number_of_snps","Frequency")
write.table(l_c,"03_hotspot_contain_snp_number_distribution_bar.txt",row.names=F,quote=F,sep="\t")

pdf("./03_barplot_hotspot_contain_snp_number_distribution.pdf",width = 4,height = 4)
p1 <-ggplot(data = l_c, mapping = aes(x =number_of_snps , y = Frequency)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.6)+
  p_theme+
  labs(x="Number of SNPs",y="Number of hotspot")+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black",angle=45,hjust=1))
print(p1)
dev.off()
#------------------------------filter top length
# filter_org <-filter(org,Log10_len < as.numeric(as.character(llen_c)))


# write.table(filter_org[,1:3],"../../../output/Tissue_merge/Cis_eQTL/01_hotspot_length_small_than_2290.bed",row.names=F,col.names =F,quote=F,sep="\t")


