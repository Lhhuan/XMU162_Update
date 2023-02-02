library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(caret)
library(stringr)
library(parallel)
library(multiROC)
library(gridExtra)
library(Hmisc) 
library(ggpubr)
library(cowplot)
library(ggsci)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    # legend.position ="none",
    plot.title = element_text(hjust = 0.5))

model <- c("AdaBoostClassifier","DecisionTreeClassifier","randomforest","RidgeClassifier","SGDClassifier","XGBClassifier")

read_AUC <-function(x){
    org<-read.table(paste0("03_",x,"_predict_AUC.txt"),header = T,sep = "\t") %>% as.data.frame()
    return(org)
}

tmp_auc <- lapply(model,read_AUC)
AUC <-do.call(rbind,tmp_auc)
colnames(AUC)[2]<-"Performance"
AUC <- mutate(AUC,variable ="AUC",.before=2)

read_file <-function(x){
    org<-read.table(paste0("03_",x,"_predict_ACC.txt"),header = T,sep = "\t") %>% as.data.frame()
    return(org)
}


tmp_re <- lapply(model,read_file)
re <- do.call(rbind,tmp_re)
re1 <-reshape2::melt(re,id="Model")
colnames(re1)[3] <- "Performance"
re2 <- bind_rows(re1,AUC)
re2$variable <- gsub("accuracy","Accuracy",re2$variable)
re2$variable <- gsub("F1_score","F1-score",re2$variable)


p1 <- ggplot(re2,mapping=aes(x=variable,y=Performance,fill=Model))+geom_bar(stat = 'identity',position=position_dodge(0.9))+ #
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    scale_fill_manual(values=c("#A593E0","#F6B352","#DE6449","#84B1ED","#ffda8e","#4f953b"))+
    p_theme+
    labs(x="",y="Performance")
pdf("./figure/031_model_select_performance.pdf",height=5,width=8)
print(p1)
dev.off()