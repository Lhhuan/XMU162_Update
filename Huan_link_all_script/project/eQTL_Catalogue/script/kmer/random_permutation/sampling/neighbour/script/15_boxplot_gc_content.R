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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/GC_content/")
org <-read.table("gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()

org$cluster <-as.factor(org$cluster)

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


p1 <- ggplot(org, aes(x=cluster, y=gc_content,fill =cluster)) + 
    geom_boxplot()+
    scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_brewer(palette="BuPu") #+
    theme_bw()+
    ggtitle("GC content")+
    p_theme +
    labs(x="Cluster",y="GC content")
pdf("15_gc_count_boxplot.pdf",height=4.7,width=4.5)
print(p1)
dev.off()
