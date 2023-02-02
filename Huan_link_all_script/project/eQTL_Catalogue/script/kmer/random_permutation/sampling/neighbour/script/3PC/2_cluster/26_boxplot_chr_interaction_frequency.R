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
library(data.table)


setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/HIC/")
sample_name=c("4DNFIQYQWPF5","4DNFIFLJLIS5.hic","4DNFI2TK7L2F","4DNFI18Q799K")
# sample_name=c("4DNFIFLJLIS5.hic","4DNFI2TK7L2F","4DNFI18Q799K")


for (sample in sample_name){
    hic <- paste0("/share/data0/QTLbase/huan/hic/",sample,"_5000.ginteractions_intra_chr_adjust.bed.gz")
    hotspot_hic <-paste0("hotspot_cluster_",sample,"_5000_intra_chr.bed.gz")
    system(paste0("bedtools intersect -a ",hic," -b ../GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_sorted.bed.gz -wo |gzip >",hotspot_hic))
    print(sample)
}
#====================================
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
    hic <- paste0("/share/data0/QTLbase/huan/hic/",sample,"_5000.ginteractions_intra_chr_adjust.bed.gz")
    hotspot_hic <-paste0("hotspot_cluster_",sample,"_5000_intra_chr.bed.gz")
    org <- fread(hic,header = F,sep = "\t") %>% as.data.frame()
    colnames(org) <-c("C_CHR","C_start","C_end","signal","line_num")
    hic_sum =org%>%group_by(C_CHR,C_start,C_end)%>%summarise(sum=n())%>%data.frame()
    mean_hic_freq <-mean(hic_sum$sum)/5
    org1 <-fread(hotspot_hic,header = F,sep = "\t") %>% as.data.frame()
    colnames(org1) <-c("C_CHR","C_start","C_end","signal","line_num","H_CHR","H_start","H_end","Cluster","overlap_bp")
    #=================
    sum_freq <-org1%>%group_by(H_CHR,H_start,H_end,Cluster)%>%summarise(sum_freq=n())%>%data.frame()
    colnames(sum_freq) <-c("CHR","start","end","cluster","sum_freq")
    # sum_freq$cluster <-as.factor(sum_freq$cluster)
    fdat <-left_join(all_hotspot,sum_freq,by=c("CHR","start","end","cluster"))
    fdat$sum_f[is.na(fdat$sum_f)] <-0
    fdat$length <- fdat$end -fdat$start
    fdat$freq_by_hotspot_length <- fdat$sum_freq/fdat$length*1000
    #==========
    fdat$cluster <-factor(fdat$cluster,levels=c(1,2))
    title = sapply(sample,function(x){link1$new[match(x,link1$old)]})
    names(title)=NULL
    type = sapply(sample,function(x){link2$new[match(x,link2$old)]})
    names(type)=NULL
    title2 <-gsub("\\s","_",title)
    # p1 <- ggplot(fdat, aes(x=cluster, y=mean_by_hotspot_length,fill =cluster)) + 
    p1 <- ggplot(fdat, aes(x=cluster, y=log(freq_by_hotspot_length),fill =cluster)) + 
        geom_boxplot()+
        scale_fill_manual(values=c("#1E77B4","#FF7F0E"))+
        theme_bw()+
        # ggtitle(paste0(title," (intra chr)"))+
        ggtitle(title)+
        geom_hline(aes(yintercept=log(mean_hic_freq)),linetype="dashed")+
        stat_compare_means(method = "wilcox.test",comparisons =list(c(1,2)),label="p.signif",method.args = list(alternative = "less"))+
        p_theme +
        labs(x="Cluster",y=paste0("Log(Frequency of interactions per kb)"))
    pdf(paste0("25_boxplot_",title2,"_interaction_frequency_by_hotspot_length_intra_chr.pdf"),height=3.1,width=3)
    print(p1)
    dev.off()
    #=====
    print(sample)
}
