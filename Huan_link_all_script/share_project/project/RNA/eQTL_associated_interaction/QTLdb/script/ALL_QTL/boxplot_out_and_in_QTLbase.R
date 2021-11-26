library(ggplot2)
library(Rcpp)
# library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gcookbook)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/")


org1<-read.table("../output/ALL_QTL/QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_QTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org2<-read.table("../output/ALL_QTL/1kg_out_QTLbase_NHPoisson_emplambda_interval_1000_cutoff_7.3_all_QTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org1$Class <-"YES"
org2$Class <- "NO"

rs<- bind_rows(org1,org2)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/ALL_QTL/")
pdf("lambda_comparison.pdf",height = 4,width = 5) 
p <- ggplot(rs,aes(x=Class,y=emplambda))+geom_boxplot(aes(fill=Class),outlier.colour = NA)+
  scale_y_continuous(name = "lambda") + scale_x_discrete(name = "Class")

p
dev.off()
print("aaa")

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/ALL_QTL/")
pdf("lambda_comparison_outlier.pdf",height = 4,width = 5) 
p <- ggplot(rs,aes(x=Class,y=emplambda))+geom_boxplot(aes(fill=Class))+
  scale_y_continuous(name = "lambda") + scale_x_discrete(name = "Class")

p
dev.off()
print("aaa")

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/ALL_QTL/")
pdf("lambda_comparison_boxplot_violin.pdf",height = 4,width = 5) 
p <- ggplot(rs,aes(x=Class,y=emplambda,fill=Class))+geom_violin(width=0.8)+geom_boxplot(width=0.5)+
  scale_y_continuous(name = "lambda") + scale_x_discrete(name = "Class")

p
dev.off()
print("aaa")


p <- ggplot(rs,aes(x=Class,y=emplambda))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ #按照标签画图
  scale_y_continuous(name = "drug sensitivity in paper") +
  scale_x_discrete(name = "top proportion")
my_comparisons <- list(c("0_0.1","0.1_0.9"), c("0_0.1","0.9_1"),c("0.1_0.9","0.9_1") ) 
p+coord_cartesian(ylim = c(-0.5,0.9))+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", method.args = list(alternative = "less"),
                                                           tip.length = 0.01,label.y=c(0.65,0.75,0.88)) ##单侧检验,假设左边小于右边 用less

##组内进行比较，指定p bar的长度并指定坐标轴p值标的位置
dev.off()