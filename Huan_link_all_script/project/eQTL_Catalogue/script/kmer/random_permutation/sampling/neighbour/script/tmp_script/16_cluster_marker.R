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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/markers/")
Cluster <-read.table("../../11_chr1_louvain_pca5_k500.txt",header = T,sep = "\t") %>% as.data.frame()
org <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/enrichment/hotspot_cutoff_0.176_marker_jaccard_index_hotspot.txt.gz",header = T,sep = "\t") %>% as.data.frame()



chr1_org <- filter(org,hotspot_chr=="chr1")
chr1_org$hotspot <-paste0(chr1_org$hotspot_chr,":",chr1_org$hotspot_start,"-",chr1_org$hotspot_end)
chr1_org<-left_join(chr1_org,Cluster[,c("cluster","hotspot")],by="hotspot")


for(marker in unique(chr1_org$Marker)){
    dat <-filter(chr1_org,Marker==marker)
    p1 <- ggplot(dat, aes(x=as.factor(cluster), y=jaacard_index,group=cluster,fill = as.factor(cluster))) + 
    geom_boxplot(outlier.colour = NA)+
    theme_bw()
    pdf(paste0("16_",marker,"_jaacard_boxplot_cluster.pdf"),height=6,width=7)
    print(p1)
    dev.off()
    print(marker)
}

chr1_org_n <- filter(chr1_org,jaacard_index>0)
chr1_orgn_marker <- chr1_org_n %>% group_by(cluster,Marker)%>%summarise(hit_number=n())%>%data.frame()
cluster_n <- Cluster%>%group_by(cluster)%>%summarise(hotspot_num=n())%>%data.frame()
chr1_orgn_marker <-left_join(chr1_orgn_marker,cluster_n,by="cluster")
chr1_orgn_marker$ratio <- chr1_orgn_marker$hit_number/chr1_orgn_marker$hotspot_num
chr1_orgn_marker$cluster <-as.factor(chr1_orgn_marker$cluster)


p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))



for(marker in unique(chr1_orgn_marker$Marker)){
    dat <-filter(chr1_orgn_marker,Marker==marker)
    marker1 <-gsub("CHROMATIN_Accessibility","CA",marker)
    p1 <-ggplot(data = dat, mapping = aes(x =cluster , y = ratio,fill=as.factor(cluster))) + geom_bar(stat = 'identity', width=0.6)+
    p_theme+scale_y_continuous(expand=c(0,0))+
    ggtitle(marker1)+
    theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black",hjust=1),
        plot.title = element_text(hjust = 0.5))+
    labs(x="Cluster",y="Fraction of hotspot")

    pdf(paste0("16_",marker,"_barplot_cluster.pdf"),height=6,width=7)
    print(p1)
    dev.off()
    print(marker)
}

#================================