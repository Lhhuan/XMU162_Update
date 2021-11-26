library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

pwd_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script"
setwd(pwd_dir)
org<-read.table("../output/merge_QTL_all_QTLtype_pop_cistrans.txt.gz",header = T,sep = "\t") %>% as.data.frame()
n_all <-nrow(org)
n_all
# 293623477
cis_1M<-filter(org,cis_trans_1MB=="cis")
n_cis_1M <-nrow(cis_1M)
n_trans_1M =n_all-n_cis_1M
n_cis_1M
# 271386276
n_trans_1M
#22237201
cis_10M<-filter(org,cis_trans_10MB=="cis")
n_cis_10M <-nrow(cis_10M)
n_trans_10M =n_all-n_cis_10M
n_cis_10M
#275994248
n_cis_10M
#275994248
n_trans_10M
# 17629229
Number<-c(n_cis_1M,n_cis_10M,n_trans_1M,n_trans_10M)
Cis_trans<-c("Cis","Cis","Trans","Trans")
Distance<-c("1M","10M","1M","10M")
a = data.frame(Number = number,Cis_trans=cis_trans,Distance =distance)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")
pdf("cis_trans_number_1M_10M.pdf",height = 5,width = 5)
p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------
p<-ggplot(a,aes(x=Distance,y=Number,fill=Cis_trans))+geom_bar(stat="identity",width=0.5)
p1<-p+geom_text(aes(label=Number,vjust=1.5))
p1<-p1+labs(fill="")
P1<-p1+p.theme
print(p1)
dev.off()