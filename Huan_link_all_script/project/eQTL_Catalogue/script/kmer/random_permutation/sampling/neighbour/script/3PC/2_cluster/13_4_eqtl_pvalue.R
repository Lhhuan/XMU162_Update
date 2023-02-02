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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/gene/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05_pvalue.bed.gz",header = F,sep = "\t") %>% as.data.frame()
Cluster <-read.table("../../11_whole_genome_leiden_pca3_k50_resolution2e-05.txt",header = T,sep = "\t") %>% as.data.frame()
all_gene$hotspot <-paste0(all_gene$V6,":",all_gene$V7,"-",all_gene$V8)
all_gene$eqtl <- paste0(all_gene$V1,":",all_gene$V2,"-",all_gene$V3)
all_gene$length <- all_gene$V8 - all_gene$V7
all_gene<-left_join(all_gene,Cluster[,c("cluster","hotspot")],by="hotspot")
colnames(all_gene)[4:5]<-c("ENSEMBL","Pvalue")
all_gene<-all_gene[,c("cluster","hotspot","eqtl","length","ENSEMBL","Pvalue")]


p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

all_gene1 <-filter(all_gene,Pvalue<0.05)

all_gene1$cluster <-factor(all_gene1$cluster,levels=c(1,2))
# p1 <- ggplot(all_gene1, aes(x=as.factor(cluster), y=-log10(Pvalue),group=cluster,fill = as.factor(cluster))) + 
p1 <- ggplot(all_gene1, aes(x=as.factor(cluster), y=Pvalue,group=cluster,fill = as.factor(cluster))) + 
    geom_boxplot()+
    # labs(x="Cluster",y="-Log10(Pvalue)")+
    labs(x="Cluster",y="p-value")+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    # scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
    # scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    scale_fill_manual(values=c("#1E77B4","#FF7F0E"))+
    theme_bw()+
    ggtitle("eQTL p-value")+
    theme(legend.position="none")+
    stat_compare_means(method = "wilcox.test",comparisons =list(c(1,2)),label="p.signif",method.args = list(alternative = "less"))+
    p_theme
pdf("13_4_eqtl_pvalue_boxplot.pdf",height=3.1,width=3)
print(p1)
dev.off()

