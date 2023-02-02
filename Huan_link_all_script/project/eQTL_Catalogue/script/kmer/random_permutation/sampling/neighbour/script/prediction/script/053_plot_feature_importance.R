setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
dat <- read.csv("05_XGB2_model_feature_importance_weight.txt",na.strings = "",sep="\t")
weight <-dat[order(-dat$Feature_importance),]
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 10),
                                                  axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))
dat$Feature <-capitalize(dat$Feature)
dat$Feature <-as.factor(dat$Feature)
dat$Feature <-gsub("_"," ",dat$Feature)
dat$Feature <-gsub("CHROMATIN Accessibility","CA",dat$Feature)
 

# pdf("./figure/10_2_not_fill_feature_importance_weight.pdf",width = 4,height = 4)
p1 <-ggplot(data = dat[1:20,], mapping = aes(x = reorder(Feature, Feature_importance), y = Feature_importance)) + geom_bar(stat = 'identity', fill = "#1A5599", width=0.8)+
  coord_flip()+p_theme+labs(y="Feature importance",x=" ",title="Weight") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))
pdf("./figure/053_top20_feature_importance_weight.pdf",width = 4,height = 4)
print(p1)
dev.off()
#
#--------------------gain
dat <- read.csv("05_XGB2_model_feature_importance_gain.txt",na.strings = "",sep="\t")
dat$Feature <-as.factor(dat$Feature)
dat$Feature <-gsub("_"," ",dat$Feature)
dat$Feature <-gsub("CHROMATIN Accessibility","CA",dat$Feature)
# pdf("./figure/10_2_not_fill_feature_importance_gain.pdf",width = 4,height = 4)
p1 <-ggplot(data = dat[1:20,], mapping = aes(x = reorder(Feature, Feature_importance), y = Feature_importance)) + geom_bar(stat = 'identity', fill = "#1A5599", width=0.8)+
  coord_flip()+p_theme+labs(y="Feature importance",x=" ",title="Gain") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))
pdf("./figure/053_top20_feature_importance_gain.pdf",width = 4,height = 4)
print(p1)
dev.off()
#---------------------shap
dat <- data.table::fread("05_XGB2_model_feature_importance_shap.txt.gz",header = T,sep = "\t")
dat <-abs(dat)
col_mean = apply(dat,2,mean)%>%as.data.frame()
col_mean$Feature <-rownames(col_mean)

colnames(col_mean) <-c("value","Feature")
# data <-sort(col_mean,descresing=T)
data1 <-col_mean[order(-col_mean$value),]
shap <-data1
data1$Feature <-capitalize(data1$Feature)
data1$Feature <-as.factor(data1$Feature)
data1$Feature <-gsub("_"," ",data1$Feature)
data1$Feature <-gsub("CHROMATIN Accessibility","CA",data1$Feature)
# data1 <-head(data,39)
# data1[40,] <- c(sum(data$value[40:539]),"Sum_of_500_other_features")

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 10),
                                                  axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))


# pdf("./figure/10_2_not_fill_feature_importance_shap.pdf",width = 4,height = 4)
p1 <-ggplot(data = data1[1:20,], mapping = aes(x =reorder(Feature,value), y = value)) + geom_bar(stat = 'identity', fill = "#1A5599", width=0.8)+
  coord_flip()+p_theme+labs(y="mean(|SHAP value|)",x=" ",title="Shap") + scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black"))
pdf("./figure/053_top20_feature_importance_shap.pdf",width = 4,height = 4)
print(p1)
dev.off()