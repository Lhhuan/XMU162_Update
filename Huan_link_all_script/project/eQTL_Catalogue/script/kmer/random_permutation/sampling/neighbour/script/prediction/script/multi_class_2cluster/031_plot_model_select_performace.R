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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/multi_class_2cluster/output/")
p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 9,angle = 45,hjust=1),#vjust=1,hjust=0.5
    axis.line = element_line(colour = "black"),
    # legend.position ="none",
    plot.title = element_text(hjust = 0.5))

model <- c("DecisionTreeClassifier","RandomForest","RidgeClassifier","XGBClassifier")


read_file <-function(x){
    org<-read.table(paste0("03_",x,"_predict_performance.txt"),header = T,sep = "\t") %>% as.data.frame()
    return(org)
}


tmp_re <- lapply(model,read_file)
re <- do.call(rbind,tmp_re)
re <-re%>%select(-c(F1_score_micro,Recall_micro,Precision_micro))
re1 <-reshape2::melt(re,id="Model")
colnames(re1)[3] <- "Performance"
# re2 <- bind_rows(re1,AUC)
re1$variable <- gsub("accuracy","Accuracy",re1$variable)
re1$variable <- gsub("F1_score","F1-score",re1$variable)
re1$variable <- gsub("_","-",re1$variable)
re1$Model <- gsub("RandomForest","RandomForestClassifier",re1$Model)

p1 <- ggplot(re1,mapping=aes(x=variable,y=Performance,fill=Model))+geom_bar(stat = 'identity',position=position_dodge(0.9))+ #
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    scale_fill_manual(values=c("#A593E0","#F6B352","#DE6449","#84B1ED"))+
    # scale_fill_manual(values=c("#A593E0","#F6B352","#DE6449","#84B1ED","#ffda8e","#4f953b"))+
    p_theme+
    labs(x="",y="Performance")
pdf("031_model_select_performance.pdf",height=5,width=8)
print(p1)
dev.off()