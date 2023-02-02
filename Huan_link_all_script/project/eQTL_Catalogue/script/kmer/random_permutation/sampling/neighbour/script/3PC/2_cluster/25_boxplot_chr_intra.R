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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/HIC/")

sample_name=c("4DNFIQYQWPF5","4DNFIFLJLIS5.hic","4DNFI2TK7L2F","4DNFI18Q799K")
# sample_name=c("4DNFIFLJLIS5.hic","4DNFI2TK7L2F","4DNFI18Q799K")


for (sample in sample_name){
    # org <-read.table(paste0("/share/data0/QTLbase/huan/hic/",sample,"_5000.ginteractions_intra_chr_adjust.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    # colnames(org) <-c("C_CHR","C_start","C_end","signal","line_number")
    # segment_score <-org%>%group_by(C_CHR,C_start,C_end)%>%summarise(segment_signal=sum(signal))%>%data.frame()
    # output_name<-paste0("/share/data0/QTLbase/huan/hic/",sample,"_5000.ginteractions_intra_chr_adjust_sum.bed")
    # write.table(segment_score,output_name,row.names = F, col.names = F,quote =F,sep="\t")
    sort_output_name <- paste0("/share/data0/QTLbase/huan/hic/",sample,"_5000.ginteractions_intra_chr_adjust_sum_sort.bed.gz")
    hotspot_hic <-paste0("hotspot_cluster_",sample,"_5000_intra_chr_adjust.bed.gz")
    # system(paste0("less ",output_name,"|sort -k1,1 -k2,2n |gzip >",sort_output_name))
    # system(paste0("gzip ",output_name))
    system(paste0("bedtools intersect -a ",sort_output_name," -b ../GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_sorted.bed.gz -wo |gzip >",hotspot_hic))
    print(sample)
}

link1 <- data.frame(old=sample_name,new=c("in situ Hi-C on H1-hESC","in situ Hi-C on HFFc6","Micro-C XL on H1-hESC","Micro-C XL on HFFc6"))
title = sapply(sample_name,function(x){link1$new[match(x,link1$old)]})
link2 <- data.frame(old=sample_name,new=c("in situ Hi-C","in situ Hi-C","Micro-C XL","Micro-C XL"))
type = sapply(sample_name,function(x){link2$new[match(x,link2$old)]})

all_hotspot <- read.table("../GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
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

for (sample in sample_name){
    sort_output_name <- paste0("/share/data0/QTLbase/huan/hic/",sample,"_5000.ginteractions_intra_chr_adjust_sum_sort.bed.gz")
    hotspot_hic <-paste0("hotspot_cluster_",sample,"_5000_intra_chr_adjust.bed.gz")
    org <- read.table(sort_output_name,header = F,sep = "\t") %>% as.data.frame()
    colnames(org) <-c("C_CHR","C_start","C_end","signal")
    mean_sig =mean(org$signal)/5
    org1 <-read.table(hotspot_hic,header = F,sep = "\t") %>% as.data.frame()
    colnames(org1) <-c("C_CHR","C_start","C_end","signal","H_CHR","H_start","H_end","Cluster","overlap_bp")
    #=================
    sum_score <-org1%>%group_by(H_CHR,H_start,H_end,Cluster)%>%summarise(sum=sum(signal))%>%data.frame()
    colnames(sum_score) <-c("CHR","start","end","cluster","sum")
    # sum_score$cluster <-as.factor(sum_score$cluster)
    fdat <-left_join(all_hotspot,sum_score,by=c("CHR","start","end","cluster"))
    fdat$sum[is.na(fdat$sum)] <-0
    fdat$length <- fdat$end -fdat$start
    fdat$mean_by_hotspot_length <- fdat$sum/fdat$length*1000
    #==========
    fdat$cluster <-factor(fdat$cluster,levels=c(1,2))
    title = sapply(sample,function(x){link1$new[match(x,link1$old)]})
    names(title)=NULL
    type = sapply(sample,function(x){link2$new[match(x,link2$old)]})
    names(type)=NULL
    title2 <-gsub("\\s","_",title)
    # p1 <- ggplot(fdat, aes(x=cluster, y=mean_by_hotspot_length,fill =cluster)) + 
    p1 <- ggplot(fdat, aes(x=cluster, y=log(mean_by_hotspot_length),fill =cluster)) + 
        geom_boxplot()+
        # scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
        scale_fill_manual(values=c("#1E77B4","#FF7F0E"))+
        theme_bw()+
        # ggtitle(paste0(title," (intra chr)"))+
        ggtitle(title)+
        geom_hline(aes(yintercept=log(mean_sig)),linetype="dashed")+
        stat_compare_means(method = "wilcox.test",comparisons =list(c(1,2)),label="p.signif",method.args = list(alternative = "less"))+
        p_theme +
        labs(x="Cluster",y=paste0("Log(mean(signal per kb))"))
    pdf(paste0("25_boxplot_",title2,"_signal_mean_by_hotspot_length_intra_chr.pdf"),height=3.1,width=3)
    print(p1)
    dev.off()
    #=====
    colnames(fdat) <-c("CHR","start","end","cluster","sum_signal","hotspot_length","mean_signal")
    org2 <- org1[,c("C_CHR","C_start","C_end","H_CHR","H_start","H_end","Cluster","overlap_bp")]%>%unique()
    sum_length <-org2 %>%group_by(H_CHR,H_start,H_end,Cluster)%>%summarise(sum_overlap=sum(overlap_bp))%>%data.frame()
    colnames(sum_length)[1:4] <-c("CHR","start","end","cluster")
    sum_length$cluster <-as.factor(sum_length$cluster)
    fdat <-left_join(fdat,sum_length,by=c("CHR","start","end","cluster"))
    fdat$sum_overlap[is.na(fdat$sum_overlap)] <-0
    fdat$mean_overlap <-fdat$sum_overlap/fdat$hotspot_length
    #===========
    p1 <- ggplot(fdat, aes(x=cluster, y=mean_overlap,fill =cluster)) + 
        geom_boxplot()+
        scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
        theme_bw()+
        # ggtitle(paste0(title," (intra chr)"))+
        ggtitle(title)+
        p_theme +
        labs(x="Cluster",y=paste0("fraction of ",type))
    pdf(paste0("25_boxplot_",title2,"_mean_overlap_by_hotspot_length_intra_chr.pdf"),height=4.7,width=4.5)
    print(p1)
    dev.off()
    print(sample)
}