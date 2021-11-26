library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)

# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase")
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")

#---------------------------------------------------------------
org<-read.table("../output/all_qtl_clump_locus.txt.gz",header = T,sep = "\t") %>% as.data.frame()

pdf("locus_region_distance_density.pdf",height = 4,width = 5) #把图片存下来
p1 <-ggplot (org,aes(x=Locus_distance,colour = QTL_type)) +geom_density(alpha=.2) + xlab("Distance of locus region") #alpha设置透明度
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Density")
# p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p1<-p1+labs(colour = "QTL type") #修改图例名字
p1
dev.off()




org<-read.table("../output/unique_qtl_locus.txt.gz",header = T,sep = "\t") %>% as.data.frame()

pdf("unique_locus_region_distance_density_log.pdf",height = 4,width = 5) #把图片存下来
p1 <-ggplot (org,aes(x=log10(Locus_distance),colour = QTL_type)) +geom_density(alpha=.2) + xlab("Length of locus region") #alpha设置透明度
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Density")
# p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p1<-p1+labs(colour = "QTL type") #修改图例名字
p1
dev.off()
# pdf("Absolute_distance_between_snp_and_Trait_density.pdf",height = 4,width = 5) #把图片存下来
# p1 <-ggplot (org,aes(x=abs_SNP.Trait_distance,colour = QTL_type)) +geom_density(alpha=.2) + xlab("Absolute distance between snp and Trait") #alpha设置透明度
# p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                                                 panel.background = element_blank(), axis.title.y = element_text(size = 10),
#                                                 axis.title.x = element_text(size = 10),
#                                                 axis.line = element_line(colour = "black")) #去背景
# p1<-p1+ylab("Density")
# # p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
# p1<-p1+labs(colour = "QTL type") #修改图例名字
# p1
# dev.off()
#----------------------------------------------geom point
png("Distance_between_snp_and_Trait_point.png",height = 4,width = 5) #把图片存下来
p2 <-ggplot (org,aes(x=SNP,y=Locus_distance, colour = QTL_type)) +geom_point(size = 0.1) + xlab("Distance between snp and Trait") 
p2<-p2+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p2<-p2+ylab("Distance")
# p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p2<-p2+labs(colour = "QTL type") #修改图例名字
p2
dev.off()
