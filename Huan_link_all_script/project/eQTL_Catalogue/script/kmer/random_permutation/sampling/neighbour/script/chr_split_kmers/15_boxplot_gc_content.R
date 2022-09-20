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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/GC_content/")
org <-read.table("gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()

org$cluster <-as.factor(org$cluster)

p1 <- ggplot(org, aes(x=cluster, y=gc_content,fill =cluster)) + 
    geom_boxplot()+
    
    # scale_fill_brewer(palette="BuPu") #+
    theme_bw()+
    theme(legend.position="none") +
pdf("15_gc_count_boxplot.pdf")
print(p1)
dev.off()


org$cluster2 <- "Other"
org[which(org$cluster==3),"cluster2"]=3
org[which(org$cluster==5),"cluster2"]=5


my_comparisons <- list( c("3", "5"), c("3", "Other"), c("5", "Other") )

p1 <- ggplot(org, aes(x=cluster2, y=gc_content,group=cluster2,fill = cluster2)) + 
    geom_boxplot()+
    stat_compare_means(comparisons = my_comparisons)+
    theme_bw()
pdf("15_gc_count_boxplot_cluster2.pdf")
print(p1)
dev.off()

org$cluster3 <- "Other"
org[which(org$cluster==3),"cluster3"]="3_5"
org[which(org$cluster==5),"cluster3"]="3_5"

p1 <- ggplot(org, aes(x=cluster3, y=gc_content,group=cluster3,fill = cluster3)) + 
    geom_boxplot()+
    stat_compare_means()+
    theme_bw()
pdf("15_gc_count_boxplot_cluster3.pdf")
print(p1)
dev.off()


org$cluster4 <- "Other"
org[which(org$cluster==3),"cluster4"]="3"
# org[which(org$cluster==5),"cluster4"]="3_5"

p1 <- ggplot(org, aes(x=cluster4, y=gc_content,group=cluster4,fill = cluster4)) + 
    geom_boxplot()+
    stat_compare_means()+
    theme_bw()
pdf("15_gc_count_boxplot_cluster4.pdf")
print(p1)
dev.off()

org$cluster5 <- "Other"
org[which(org$cluster==5),"cluster5"]="5"
# org[which(org$cluster==5),"cluster4"]="3_5"

p1 <- ggplot(org, aes(x=cluster5, y=gc_content,group=cluster5,fill = cluster5)) + 
    geom_boxplot()+
    stat_compare_means()+
    theme_bw()
pdf("15_gc_count_boxplot_cluster5.pdf")
print(p1)
dev.off()

org$cluster6 <- "Other"
org[which(org$cluster==2),"cluster6"]="2"

p1 <- ggplot(org, aes(x=cluster6, y=gc_content,group=cluster6,fill = cluster6)) + 
    geom_boxplot()+
    stat_compare_means()+
    theme_bw()
pdf("15_gc_count_boxplot_cluster6.pdf")
print(p1)
dev.off()