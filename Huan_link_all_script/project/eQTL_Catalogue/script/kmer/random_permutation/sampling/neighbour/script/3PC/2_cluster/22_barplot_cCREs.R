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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/Cis_Regulatory_Elements/")
org <-read.table("hotspot_cluster_cCREs.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)<-c("CHR","start","end","cluster","cre_chr","cre_start","cre_end","name","score","strand","thickStart","thickEnd","reserved","ccre","encodeLabel","zScore","ucscLabel","accessionLabel","description")
org$hotspot <-paste(org$CHR,org$start,org$end,sep=":")
org$cluster <-as.factor(org$cluster)
org <- org[,c("hotspot","cluster","ucscLabel","name")]
org1 <- unique(org[,c("hotspot","cluster","ucscLabel")])

dat <-org1%>%group_by(ucscLabel,cluster)%>%summarise(ucscLabel_count=n())%>%data.frame()

all_hotspot <- read.table("../GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(all_hotspot) <-c("CHR","start","end","cluster")
all_hotspot$cluster <-as.factor(all_hotspot$cluster)
all_hotspot$hotspot <- paste(all_hotspot$CHR,all_hotspot$start,all_hotspot$end,sep=":")
all_hotspot <- all_hotspot[,c("hotspot","cluster")]
cluster_n <- all_hotspot %>%group_by(cluster)%>%summarise(hotspot_count=n())%>%data.frame()


dat <- left_join(dat,cluster_n,by=c("cluster"))
dat$proporation <-dat$ucscLabel_count/dat$hotspot_count

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

plist<-list()
for(i in unique(dat$ucscLabel)){
    fdatcount <-filter(dat,ucscLabel==i)
    fdatcount$cluster <-factor(fdatcount$cluster,levels=c(1,2))
   plist[[i]] <- ggplot(fdatcount, aes(x=as.factor(cluster), y=proporation,group=cluster,fill = as.factor(cluster))) + 
    # geom_boxplot(outlier.colour = NA)+
    geom_bar(stat = 'identity', width=0.8)+
    # scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    scale_fill_manual(values=c("#1E77B4","#FF7F0E"))+
    theme_bw()+
    labs(x="Cluster",y="Fraction of hotspots")+
    ggtitle(i)+
    theme(legend.position ="none")+
    p_theme
    print(i)
    pdf(paste0("12_barplot_of_",i,".pdf"),height=4.7,width=4.5)
    print(plist[[i]])
    dev.off()
}

pdf("22_barplot_of_cCRE.pdf",width=11, height=2.8)
gridExtra::marrangeGrob(plist,ncol=5,nrow=1,top="Candidate Cis-Regulatory Elements")
dev.off()


