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

feature <- c("GC_content","Kmer") #,"marker"

read_AUC <-function(x){
    org<-read.table(paste0("051_",x,"_XGBClassifier_final_parameter_predict_AUC.txt"),header = T,sep = "\t") %>% as.data.frame()
    return(org)
}

tmp_auc <- lapply(feature,read_AUC)
AUC <-do.call(rbind,tmp_auc)
ALL_AUC <- read.table("05_XGBClassifier_final_parameter_predict_AUC.txt",header = T,sep = "\t") %>% as.data.frame()
ALL_AUC[1,1] <-"ALL features"
colnames(ALL_AUC)[1] <- "Feature_type"
AUC <-bind_rows(AUC,ALL_AUC)
colnames(AUC)[2]<-"Performance"
AUC <- mutate(AUC,variable ="AUC",.before=2)

read_file <-function(x){
    org<-read.table(paste0("051_",x,"_XGBClassifier_final_parameter_predict_ACC.txt"),header = T,sep = "\t") %>% as.data.frame()
    return(org)
}


tmp_re <- lapply(feature,read_file)
re <- do.call(rbind,tmp_re)
ALL_ACC <-read.table("05_XGBClassifier_final_parameter_predict_ACC_without_CTCF.txt",header = T,sep = "\t") %>% as.data.frame()
ALL_ACC[1,1] <-"ALL features"
colnames(ALL_ACC)[1] <- "Feature_type"
re <- bind_rows(re,ALL_ACC)
#=============
# marker <- read.table("051_marker_XGBClassifier_final_parameter_predict_ACC_without_ctcf.txt",header = T,sep = "\t") %>% as.data.frame()
# ALL_ACC[1,1] <-"ALL features"
# colnames(ALL_ACC)[1] <- "Feature_type"
# re <- bind_rows(re,ALL_ACC)
# #===================
re1 <-reshape2::melt(re,id="Feature_type")
colnames(re1)[3] <- "Performance"
re2 <- bind_rows(re1,AUC)
re2$variable <- gsub("accuracy","Accuracy",re2$variable)
re2$variable <- gsub("F1_score","F1-score",re2$variable)
re2$Feature_type <- gsub("Markers","Marker",re2$Feature_type)
re2$Feature_type <- gsub("Marker","Markers",re2$Feature_type)
re2$Feature_type <- gsub("GC_content","GC content",re2$Feature_type)
p1 <- ggplot(re2,mapping=aes(x=variable,y=Performance,fill=Feature_type))+geom_bar(stat = 'identity',position=position_dodge(0.9))+ #
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    scale_fill_manual(values=c("#A593E0","#F6B352","#DE6449","#84B1ED"))+
    p_theme+
    labs(x="",y="Performance",fill="Feature type")
pdf("./figure/052_feature_select_performance.pdf",height=5,width=6.5)
print(p1)
dev.off()