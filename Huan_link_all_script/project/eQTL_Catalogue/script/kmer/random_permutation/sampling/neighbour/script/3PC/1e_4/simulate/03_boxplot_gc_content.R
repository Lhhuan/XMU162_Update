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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/3PC/1e_4/simulate/output/GC_content/")
C1 <-read.table("01_random_c1_gc_content.bed.gz",header = F,sep = "\t") %>% as.data.frame()
C2 <-read.table("01_random_c2_gc_content.bed.gz",header = F,sep = "\t") %>% as.data.frame()
C3 <-read.table("01_random_c3_gc_content.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(C1)[5] <- "gc_content"
colnames(C2)[5] <- "gc_content"
colnames(C3)[5] <- "gc_content"
C1$cluster <-"1"
C2$cluster <-"2"
C3$cluster <-"3"
org <- bind_rows(C1,C2,C3)
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
# org$cluster <-factor(org$cluster,levels=c(3,1,2))
p1 <- ggplot(org, aes(x=cluster, y=gc_content,fill =cluster)) + 
    geom_boxplot()+
    # scale_fill_manual(values=c("#2CA02C","#1E77B4","#FF7F0E"))+
    scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C"))+
    # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
    # scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
    # scale_fill_manual(breaks=c(6,4,2,1,5,3),values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    # scale_fill_brewer(palette="BuPu") #+
    theme_bw()+
    ggtitle("GC content")+
    p_theme +
    labs(x="Cluster",y="GC content")
pdf("03_gc_count_boxplot.pdf",height=4.7,width=4.5)
print(p1)
dev.off()

