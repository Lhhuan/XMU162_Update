library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(parallel)



# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/")
# org<-read.table("09_prepare_TPR_FPR_ROC_chr1_22.txt",header = T,sep = "\t") %>% as.data.frame()
# enhancer <-dplyr::filter(org,Factor == "enhancer")
# promoter <-dplyr::filter(org,Factor == "promoter")
# pdf("enhancer.pdf")
# p <-ggplot(enhancer,mapping = aes(x = FPR, y = TPR,colour = cutoff)) + geom_point(size = 3)
# print(p)
# dev.off()

# pdf("promoter.pdf")
# p <-ggplot(promoter,mapping = aes(x = FPR, y = TPR,colour = cutoff)) + geom_point(size = 3)
# print(p)
# dev.off()




setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/ALL/")
org<-read.table("08_prepare_number_ROC_refine_factor_count.txt",header = T,sep = "\t") %>% as.data.frame()
enhancer <-dplyr::filter(org,Factor == "enhancer")
promoter <-dplyr::filter(org,Factor == "promoter")
pdf("enhancer_factor_refine.pdf")
p <-ggplot(enhancer,mapping = aes(x = FPR, y = TPR,colour = Cutoff)) + geom_point(size = 3)
print(p)
dev.off()

pdf("promoter_factor_refine.pdf")
p <-ggplot(promoter,mapping = aes(x = FPR, y = TPR,colour = Cutoff)) + geom_point(size = 3)
print(p)
dev.off()




setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/ROC/interval_18/chr1_22/")
org<-read.table("09_prepare_TPR_FPR_ROC_chr1_22.txt",header = T,sep = "\t") %>% as.data.frame()
enhancer <-dplyr::filter(org,Factor == "enhancer")
promoter <-dplyr::filter(org,Factor == "promoter")
# enhancer_tmp1 <-data.frame(cutoff =0,Factor="enhancer",TPR=0,FPR=0)
# enhancer_tmp2 <-data.frame(cutoff =1,Factor="enhancer",TPR=1,FPR=1)
# enhancer <-bind_rows(enhancer,enhancer_tmp1)
# enhancer <-bind_rows(enhancer,enhancer_tmp2)
pdf("enhancer.pdf")
p <-ggplot(enhancer,mapping = aes(x = FPR, y = TPR,colour = cutoff)) + geom_point(size = 3)
print(p)
dev.off()

pdf("promoter.pdf")
p <-ggplot(promoter,mapping = aes(x = FPR, y = TPR,colour = cutoff)) + geom_point(size = 3)
print(p)
dev.off()