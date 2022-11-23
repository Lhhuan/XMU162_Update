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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/GC_content/")
org <-read.table("gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()



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

# org$cluster <-as.factor(org$cluster)
org$cluster <-factor(org$cluster,levels=c(4,2,1,5,6,3))
p1 <- ggplot(org, aes(x=cluster, y=gc_content,fill =cluster)) + 
    geom_boxplot()+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
    # scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
    scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    # scale_fill_manual(breaks=c(6,4,2,1,5,3),values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    # scale_fill_brewer(palette="BuPu") #+
    theme_bw()+
    ggtitle("GC content")+
    p_theme +
    labs(x="Cluster",y="GC content")
pdf("15_gc_count_boxplot.pdf",height=3.1,width=3)
print(p1)
dev.off()

