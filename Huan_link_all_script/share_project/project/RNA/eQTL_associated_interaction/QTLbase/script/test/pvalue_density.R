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

org1<-read.table("../output/QTL_specific_Pvalue.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# org2<-read.table("../output/unique_by_qtl_pop_locus_r_square0.5.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# org3<-read.table("../output/unique_by_qtl_pop_locus_r_square0.8.txt.gz",header = T,sep = "\t") %>% as.data.frame()





pdf("P_vlaue_density.pdf",height = 4,width = 5) #把图片存下来
p1 <-ggplot (org1,aes(x=log10(Pvalue),colour = QTL_type)) +geom_density(alpha=.2) #alpha设置透明度
p1<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Density")+xlab("log10(Pvalue)")
p1<-p1+labs(colour = "QTL type") #修改图例名字
p1<-p1+ggtitle("Pvalue density")+theme(plot.title=element_text(size=9,hjust = 0.5)) 
p1
dev.off()
