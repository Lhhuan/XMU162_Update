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
library(ggsci)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/HIC/")

sample_name=c("4DNFIQYQWPF5","4DNFIFLJLIS5.hic","4DNFI2TK7L2F","4DNFI18Q799K")

link1 <- data.frame(old=sample_name,new=c("in situ Hi-C on H1-hESC","in situ Hi-C on HFFc6","Micro-C XL on H1-hESC","Micro-C XL on HFFc6"))
title = sapply(sample_name,function(x){link1$new[match(x,link1$old)]})
link2 <- data.frame(old=sample_name,new=c("in situ Hi-C","in situ Hi-C","Micro-C XL","Micro-C XL"))
type = sapply(sample_name,function(x){link2$new[match(x,link2$old)]})
# title2 <-gsub("\\s","_",title)
# for(sample in sample_name){
#     title <-gsub("4DNFI2TK7L2F","Micro-C XL on H1-hESC",sample)
#     title <-gsub("4DNFI18Q799K","Micro-C XL on HFFc6",sample)
#     title <-gsub("4DNFIFLJLIS5.hic","in situ Hi-C on HFFc6",sample)
#     title <-gsub("4DNFIQYQWPF5","in situ Hi-C on H1-hESC",sample)
#     type <- gsub("(4DNFIQYQWPF5)|(4DNFIFLJLIS5.hic)","in situ Hi-C",sample)
#     type <- gsub("(4DNFI2TK7L2F)|(4DNFI18Q799K)","Micro-C XL",sample)
#     title2 <-gsub("\\s","_",title)
#     print(c(sample, title,title2))
# }

for (sample in sample_name){
    org <-read.table(paste0("hotspot_cluster_",sample,"_5000_intra_chr.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    colnames(org) <-c("C_CHR","C_start","C_end","signal","H_CHR","H_start","H_end","Cluster","overlap_bp")

    all_hotspot <- read.table("../GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
    colnames(all_hotspot) <-c("CHR","start","end","cluster")

    p_theme<-theme(
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        legend.position ="none",
        plot.title = element_text(hjust = 0.5))


    sum_score <-org%>%group_by(H_CHR,H_start,H_end,Cluster)%>%summarise(sum=sum(signal))%>%data.frame()
    colnames(sum_score) <-c("CHR","start","end","cluster","sum")
    # sum_score$cluster <-as.factor(sum_score$cluster)
    all_hotspot <-left_join(all_hotspot,sum_score,by=c("CHR","start","end","cluster"))
    all_hotspot$sum[is.na(all_hotspot$sum)] <-0
    all_hotspot$length <- all_hotspot$end -all_hotspot$start
    all_hotspot$mean_by_hotspot_length <- all_hotspot$sum/all_hotspot$length

    all_hotspot$cluster <-factor(all_hotspot$cluster,levels=c(4,2,1,5,6,3))
    title = sapply(sample,function(x){link1$new[match(x,link1$old)]})
    names(title)=NULL
    type = sapply(sample,function(x){link2$new[match(x,link2$old)]})
    names(type)=NULL
    title2 <-gsub("\\s","_",title)
    # p1 <- ggplot(all_hotspot, aes(x=cluster, y=mean_by_hotspot_length,fill =cluster)) + 
    p1 <- ggplot(all_hotspot, aes(x=cluster, y=log(mean_by_hotspot_length),fill =cluster)) + 
        geom_boxplot()+
        # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
        # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
        scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
        theme_bw()+
        ggtitle(paste0(title," (intra chr)"))+
        p_theme +
        # labs(x="Cluster",y="mean(signal of hic)")
        labs(x="Cluster",y=paste0("Log(mean(signal of ",type,"))"))
    pdf(paste0("25_boxplot_",title2,"_signal_mean_by_hotspot_length_intra_chr.pdf"),height=4.7,width=4.5)
    print(p1)
    dev.off()

    colnames(all_hotspot) <-c("CHR","start","end","cluster","sum_signal","hotspot_length","mean_signal")
    org2 <- org[,c("C_CHR","C_start","C_end","H_CHR","H_start","H_end","Cluster","overlap_bp")]%>%unique()
    sum_length <-org2 %>%group_by(H_CHR,H_start,H_end,Cluster)%>%summarise(sum_overlap=sum(overlap_bp))%>%data.frame()
    colnames(sum_length)[1:4] <-c("CHR","start","end","cluster")
    sum_length$cluster <-as.factor(sum_length$cluster)
    all_hotspot <-left_join(all_hotspot,sum_length,by=c("CHR","start","end","cluster"))
    all_hotspot$sum_overlap[is.na(all_hotspot$sum_overlap)] <-0
    all_hotspot$mean_overlap <-all_hotspot$sum_overlap/all_hotspot$hotspot_length
    
    p1 <- ggplot(all_hotspot, aes(x=cluster, y=mean_overlap,fill =cluster)) + 
        geom_boxplot()+
        # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
        # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
        scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
        theme_bw()+
        ggtitle(paste0(title," (intra chr)"))+
        p_theme +
        # labs(x="Cluster",y="mean(signal of hic)")
        labs(x="Cluster",y=paste0("fraction of ",type))
    pdf(paste0("25_boxplot_",title2,"_mean_overlap_by_hotspot_length_intra_chr.pdf"),height=4.7,width=4.5)
    print(p1)
    dev.off()
}