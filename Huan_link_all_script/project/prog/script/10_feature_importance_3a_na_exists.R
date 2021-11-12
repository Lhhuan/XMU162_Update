
setwd("/home/huanhuan/project/prog/output/")
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
dat <- read.csv("06_step1_feature_importance_3a_not_fill.txt",na.strings = "",sep="\t")

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 8),
                                                  axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))

dat$feature <-as.factor(dat$feature)
dat$feature <-gsub("_"," ",dat$feature)
dat$feature <-capitalize(dat$feature) 

pdf("./figure/10_feature_importance_3a_na_exist_step1.pdf",width = 5,height = 6)
p1 <-ggplot(data = dat, mapping = aes(x = reorder(feature, weight), y = weight*100)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.5)+
  coord_flip()+p_theme+labs(y="Feature importance(%)",x=" ",title="Step1") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))
print(p1)
dev.off()
#--------------------step6
dat <- read.csv("06_step6_feature_importance_3a_not_fill.txt",na.strings = "",sep="\t")
dat$feature <-as.factor(dat$feature)
dat$feature <-gsub("_"," ",dat$feature)
dat$feature <-capitalize(dat$feature) 
pdf("./figure/10_feature_importance_3a_na_exist_step6.pdf",width = 5,height = 6)
p1 <-ggplot(data = dat, mapping = aes(x = reorder(feature, weight), y = weight*100)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.5)+
  coord_flip()+p_theme+labs(y="Feature importance(%)",x=" ",title="Step6") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))
print(p1)
dev.off()