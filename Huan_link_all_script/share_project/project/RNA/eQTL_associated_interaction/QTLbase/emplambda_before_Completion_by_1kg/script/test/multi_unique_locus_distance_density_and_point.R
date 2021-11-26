library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase")
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")

org1<-read.table("../output/unique_by_qtl_pop_locus_r_square0.3.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org2<-read.table("../output/unique_by_qtl_pop_locus_r_square0.5.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org3<-read.table("../output/unique_by_qtl_pop_locus_r_square0.8.txt.gz",header = T,sep = "\t") %>% as.data.frame()





pdf("unique_by_qtl_pop_locus_region_distance_density_log10_r_square0.3.pdf",height = 4,width = 5) #把图片存下来
p1 <-ggplot (org1,aes(x=log10(Locus_distance),colour = QTL_type)) +geom_density(alpha=.2) + xlab("log10(length of locus region)")  #alpha设置透明度
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Density")
p1<-p1+labs(colour = "QTL type") #修改图例名字
p1<-p1+ggtitle("R-square: 0.3")+theme(plot.title=element_text(size=9,hjust = 0.5)) 
p1
dev.off()

#----------------------------------------------------------
pdf("unique_by_qtl_pop_locus_region_distance_density_log10_r_square0.5.pdf",height = 4,width = 5) 
p1 <-ggplot (org2,aes(x=log10(Locus_distance),colour = QTL_type)) +geom_density(alpha=.2) + xlab("log10(length of locus region)") 
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) 
p1<-p1+ylab("Density")
p1<-p1+labs(colour = "QTL type") 
p1<-p1+ggtitle("R-square: 0.5")+theme(plot.title=element_text(size=9,hjust = 0.5)) 
p1
dev.off()

pdf("unique_by_qtl_pop_locus_region_distance_density_r_square0.5.pdf",height = 4,width = 5) 
p1 <-ggplot (org2,aes(x=Locus_distance,colour = QTL_type)) +geom_density(alpha=.2) + xlab("log10(length of locus region)") 
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) 
p1<-p1+ylab("Density")
p1<-p1+labs(colour = "QTL type") 
p1<-p1+ggtitle("R-square: 0.5")+theme(plot.title=element_text(size=9,hjust = 0.5)) 
p1
dev.off()

#-----------------------------
pdf("unique_by_qtl_pop_locus_region_distance_density_log10_r_square0.8.pdf",height = 4,width = 5) #把图片存下来
p1 <-ggplot (org3,aes(x=log10(Locus_distance),colour = QTL_type)) +geom_density(alpha=.2) + xlab("log10(length of locus region)")  #alpha设置透明度
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Density")
p1<-p1+labs(colour = "QTL type") #修改图例名字
p1<-p1+ggtitle("R-square: 0.8")+theme(plot.title=element_text(size=9,hjust = 0.5)) 
p1
dev.off()

#-------------------------------



org1$Class <-c("R-square: 0.3")
org2$Class <-c("R-square: 0.5")
org3$Class <-c("R-square: 0.8")

org<- bind_rows(org1,org2,org3)
type <-org1 %>% dplyr::select(class)%>%unique()



# #-------------------------------------------------------------------------------------------------------------------

# png("Unique_distance_between_snp_and_Trait_point.png",height = 4,width = 5) #把图片存下来
# p2 <-ggplot (org2,aes(x=Locus_region,y=Locus_distance, colour = QTL_type)) +geom_point(size = 0.1) + xlab("Distance between snp and Trait") 
# p2<-p2+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                                                 panel.background = element_blank(), axis.title.y = element_text(size = 10),
#                                                 axis.title.x = element_text(size = 10),
#                                                 axis.line = element_line(colour = "black")) #去背景
# p2<-p2+ylab("Distance")
# # p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
# p2<-p2+labs(colour = "QTL type") #修改图例名字
# p2
# dev.off()


# pdf("Unique_distance_between_snp_and_Trait_point.pdf",height = 4,width = 5) #把图片存下来
# p2 <-ggplot (org2,aes(x=Locus_region,y=Locus_distance, colour = QTL_type)) +geom_point(size = 0.1) + xlab("Distance between snp and Trait") 
# p2<-p2+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                                                 panel.background = element_blank(), axis.title.y = element_text(size = 10),
#                                                 axis.title.x = element_text(size = 10),
#                                                 axis.line = element_line(colour = "black")) #去背景
# p2<-p2+ylab("Distance")
# # p1<-p1+scale_colour_discrete(breaks= c("SNV/Indel","CNV", "Deletion","Duplication","Inversion","Translocation")) #修改图例顺序
# p2<-p2+labs(colour = "QTL type") #修改图例名字
# p2
# dev.off()