library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase")

org<-read.table("./output/count_snp_trait_pos_distance.txt",header = T,sep = "\t") %>% as.data.frame()
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")
pdf("Distance_between_snp_and_Trait_density.pdf",height = 4,width = 5) #把图片存下来
p1 <-ggplot (org,aes(x=SNP.Trait_distance,colour = QTL_type)) +geom_density(alpha=.2) + xlab("Distance between snp and Trait") #alpha设置透明度
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Density")
# p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p1<-p1+labs(colour = "QTL type") #修改图例名字
p1
dev.off()

pdf("Absolute_distance_between_snp_and_Trait_density.pdf",height = 4,width = 5) #把图片存下来
p1 <-ggplot (org,aes(x=abs_SNP.Trait_distance,colour = QTL_type)) +geom_density(alpha=.2) + xlab("Absolute distance between snp and Trait") #alpha设置透明度
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Density")
# p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p1<-p1+labs(colour = "QTL type") #修改图例名字
p1
dev.off()
#----------------------------------------------geom
pdf("Distance_between_snp_and_Trait_point.pdf",height = 4,width = 5) #把图片存下来
p2 <-ggplot (org,aes(x=SNP_chr,y=SNP.Trait_distance, colour = QTL_type)) +geom_point() + xlab("Distance between snp and Trait") 
p2<-p2+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p2<-p2+ylab("Distance")
# p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p2<-p2+labs(colour = "QTL type") #修改图例名字
p2
dev.off()


pdf("Absolute_distance_between_snp_and_Trait_geom.pdf",height = 3.5,width = 5) #把图片存下来
p1 <-ggplot (data = org, mapping = aes(x = x, y = abs_SNP.Trait_distance)) + geom_point(size = 3) + xlab("QTL") #alpha设置透明度
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 13),
                                                axis.title.x = element_text(size = 13),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Absolute distance between snp and Trait")
# p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
p1<-p1+labs(colour = "QTL type") #修改图例名字
p1
dev.off()