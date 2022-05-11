setwd("/home/huanhuan/project/prog/output/")
library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
dat <- read.csv("/home/huanhuan/project/prog/output/10_1_not_fill_feature_importance_weight.txt",na.strings = "",sep="\t")
weight <-dat[order(-dat$Feature_importance),]
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 8),
                                                  axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))

dat$Feature <-as.factor(dat$Feature)
dat$Feature <-gsub("_"," ",dat$Feature)
dat$Feature <-gsub("Lym Mono","Lym/Mono",dat$Feature)
dat$Feature <-capitalize(dat$Feature) 

pdf("./figure/10_2_not_fill_feature_importance_weight.pdf",width = 4,height = 4)
p1 <-ggplot(data = dat, mapping = aes(x = reorder(Feature, Feature_importance), y = Feature_importance)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.5)+
  coord_flip()+p_theme+labs(y="Feature importance",x=" ",title="Weight") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black",size=6),
    axis.text.x = element_text(color="black"))
print(p1)
dev.off()
weight_a <-dat
weight_a$Ratio <- weight_a$Feature_importance/sum(weight_a$Feature_importance)
weight_a$Ratio_adjust1 <- round(weight_a$Ratio *100,2)
weight_a$Ratio_adjust2 <- round(weight_a$Feature_importance /5)
write.table(weight_a,"10_2_feature_ratio.txt",col.names=T,row.names=F,quote=F,sep="\t")
#
#--------------------gain
dat <- read.csv("/home/huanhuan/project/prog/output/10_1_not_fill_feature_importance_gain.txt",na.strings = "",sep="\t")
dat$Feature <-as.factor(dat$Feature)
dat$Feature <-gsub("_"," ",dat$Feature)
dat$Feature <-gsub("Lym Mono","Lym/Mono",dat$Feature)
dat$Feature <-capitalize(dat$Feature) 

pdf("./figure/10_2_not_fill_feature_importance_gain.pdf",width = 4,height = 4)
p1 <-ggplot(data = dat, mapping = aes(x = reorder(Feature, Feature_importance), y = Feature_importance)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.5)+
  coord_flip()+p_theme+labs(y="Feature importance",x=" ",title="Gain") +scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black",size=6),
    axis.text.x = element_text(color="black"))
print(p1)
dev.off()

#---------------------shap
dat <- read.csv("/home/huanhuan/project/prog/output/10_1_not_fill_feature_importance_shap_class.txt",na.strings = "",sep="\t")
dat <-abs(dat)
col_mean = apply(dat,2,mean)%>%as.data.frame()
col_mean$Feature <-rownames(col_mean)

colnames(col_mean) <-c("value","Feature")
# data <-sort(col_mean,descresing=T)
data1 <-col_mean[order(-col_mean$value),]
shap <-data1
data1$Feature <-as.factor(data1$Feature)
data1$Feature <-gsub("_"," ",data1$Feature)
data1$Feature <-gsub("Lym Mono","Lym/Mono",data1$Feature)
data1$Feature <-capitalize(data1$Feature) 
data1$value <-as.numeric(data1$value)

# data1 <-head(data,39)
# data1[40,] <- c(sum(data$value[40:539]),"Sum_of_500_other_features")

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                  axis.title.y = element_text(size = 8),
                                                  axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))


pdf("./figure/10_2_not_fill_feature_importance_shap.pdf",width = 4,height = 4)
p1 <-ggplot(data = data1, mapping = aes(x =reorder(Feature,value), y = value)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.5)+
  coord_flip()+p_theme+labs(y="mean(|SHAP value|)",x=" ",title="Shap") + scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black",size=6),
    axis.text.x = element_text(color="black"))
print(p1)
dev.off()

#------------------
# overlap_f <- intersect(c(head(weight$Feature,15)),c(head(shap$Feature,15)))
# save(overlap_f,file="shap_weight_overlap_feature.Rdata")