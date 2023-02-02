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

setwd("/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/signalValue/hotspot/mean")
hotspot_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_max_mean_median_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz")
    org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
    org <-org[,c(1:4)]
    org$hotspot <-paste0(org$Chr,":",org$start,"-",org$end)
    cc <- read.table("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/11_whole_genome_leiden_pca3_k50_resolution1e-04.txt",header = T,sep = "\t") %>% as.data.frame()
    cc <- left_join(cc[,c("cluster","hotspot")],org[,c("hotspot","mean_signalvalue")],by="hotspot")
    cc$marker <-marker
    cc[is.na(cc)] <-0
    # return(cluster)
    marker <-gsub("CHROMATIN_Accessibility", "CA",marker)
    #=================
    cc$cluster <-factor(cc$cluster,levels=c(4,2,1,5,6,3))
    p1 <- ggplot(cc, aes(x=cluster, y=mean_signalvalue,fill =cluster)) + 
        geom_boxplot()+
        scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
        # scale_fill_manual(breaks=c(6,4,2,1,5,3),values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
        # scale_fill_brewer(palette="BuPu") #+
        theme_bw()+
        ggtitle(marker)+
        p_theme +
        labs(x="Cluster",y="Signalvalue")
    pdf(paste0("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/markers/16_",marker,"_boxplot.pdf"),height=3.1,width=3)
    print(p1)
    dev.off()
    return(p1)
}
markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","CTCF","TFBS","CHROMATIN_Accessibility")
plist <-lapply(markers,hotspot_signal)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/markers/")
pdf("16_2_marker_combine.pdf",width=8.3, height=7)
CombinePlots(plist,ncol=4,nrow=3)
dev.off()
