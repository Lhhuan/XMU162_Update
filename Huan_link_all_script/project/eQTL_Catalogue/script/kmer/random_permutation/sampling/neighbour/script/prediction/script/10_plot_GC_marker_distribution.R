library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)
library(reshape2)
library(R.utils)
library(mclust)
library(umap)
library(uwot)
library(kernlab)
library(data.table)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
org<-fread("08_warm_region_predict.txt",header = T,sep = "\t") %>% as.data.frame()

key_org <- org[,c(1:5,402:412)]
key_org$predict_class_lable <- key_org$predict_class
key_org$predict_class_lable <-gsub("1","True",key_org$predict_class_lable)
key_org$predict_class_lable <-gsub("0","False",key_org$predict_class_lable)

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))


key_org$predict_class_lable <-factor(key_org$predict_class_lable,levels=c("True","False"))
p1 <- ggplot(key_org, aes(x=as.factor(predict_class_lable), y=GC_content,group=predict_class_lable,fill = as.factor(predict_class_lable))) + 
geom_boxplot()+
scale_fill_manual(values=c("#A593E0","#F68657"))+
theme_bw()+
labs(x=NULL,y="GC content")+
ggtitle("GC content")+
theme(legend.position ="none")+
# stat_compare_means(method = "wilcox.test",label="p.signif")+
stat_compare_means(method = "wilcox.test",comparisons =list(c("True","False")),label="p.signif")+
p_theme

pdf("./figure/10_predicted_CG_content.pdf",height=2,width=2.1)
print(p1)
dev.off()


key_org1 <- key_org[,c("hotspot","predict_class_lable","CHROMATIN_Accessibility","TFBS","CTCF","H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3")]

markers_plot <-function(i){
    kk <-key_org1
    aa <-colnames(kk)[i]
    aa <-gsub("CHROMATIN_Accessibility", "CA",aa)
    colnames(kk)[i] <- "BB"
   p1 <- ggplot(kk, aes(x=as.factor(predict_class_lable), y=BB,group=predict_class_lable,fill = as.factor(predict_class_lable))) + 
    geom_boxplot()+
    scale_fill_manual(values=c("#A593E0","#F68657"))+
    theme_bw()+
    labs(x=NULL,y="Signal value")+
    ggtitle(aa)+
    theme(legend.position ="none")+
    stat_compare_means(method = "wilcox.test",comparisons =list(c("True","False")),label="p.signif")+
    p_theme

    pdf(paste0("./figure/10_predicted_",aa,".pdf"),height=3.1,width=3)
    print(p1)
    dev.off()
    print(aa)
    return(p1)
}
plist <- lapply(3:12,markers_plot)


pdf("./figure/10_marker_combine.pdf",width=8.3, height=6)
CombinePlots(plist,ncol=4,nrow=3)
dev.off()

#==============CTCF
p1 <- ggplot(key_org, aes(x=as.factor(predict_class_lable), y=GC_content,group=predict_class_lable,fill = as.factor(predict_class_lable))) + 
geom_boxplot()+
scale_fill_manual(values=c("#A593E0","#F68657"))+
theme_bw()+
labs(x=NULL,y="Signal value")+
ggtitle("CTCF")+
ylim(0,1.1)+
theme(legend.position ="none")+
# stat_compare_means(method = "wilcox.test",label="p.signif")+
stat_compare_means(method = "wilcox.test",comparisons =list(c("True","False")),label="p.signif")+
p_theme

pdf("./figure/10_predicted_CTCF_Y_limit.pdf",height=2,width=2.1)
print(p1)
dev.off()
#=============








markers_plot_log <-function(i){
    kk <-key_org
    aa <-colnames(kk)[i]
    aa <-gsub("CHROMATIN_Accessibility", "CA",aa)
    colnames(kk)[i] <- "BB"
   p1 <- ggplot(kk, aes(x=as.factor(predict_class_lable), y=log(BB+0.0000001),group=predict_class_lable,fill = as.factor(predict_class_lable))) + 
    geom_boxplot()+
    scale_fill_manual(values=c("#A593E0","#F68657"))+
    theme_bw()+
    labs(x="",y="Signal value")+
    ggtitle(aa)+
    theme(legend.position ="none")+
    stat_compare_means(method = "wilcox.test",comparisons =list(c("True","False")),label="p.signif")+
    p_theme

    pdf(paste0("./figure/10_predicted_",aa,"_log.pdf"),height=3.1,width=3)
    print(p1)
    dev.off()
    print(aa)
}
lapply(6:15,markers_plot_log )