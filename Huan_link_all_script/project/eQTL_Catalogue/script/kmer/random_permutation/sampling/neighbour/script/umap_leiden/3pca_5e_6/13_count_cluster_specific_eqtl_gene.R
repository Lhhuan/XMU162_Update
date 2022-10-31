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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/3pca_5e_6/gene/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
Cluster <-read.table("../../../11_whole_genome_umap_leiden_pca3_k30_resolution5e-06.txt",header = T,sep = "\t") %>% as.data.frame()
all_gene$hotspot <-paste0(all_gene$V1,":",all_gene$V2,"-",all_gene$V3)
all_gene$eqtl <- paste0(all_gene$V4,":",all_gene$V5,"-",all_gene$V6)
all_gene$length <- all_gene$V3 - all_gene$V2
all_gene<-left_join(all_gene,Cluster[,c("cluster","hotspot")],by="hotspot")
colnames(all_gene)[7]<-"ENSEMBL"
all_gene<-all_gene[,c("cluster","hotspot","eqtl","length","ENSEMBL")]
all_gene_count <-all_gene%>%group_by(hotspot,cluster,length)%>%summarise(eqtl_gene_count=n())%>%data.frame()
all_gene_count$adjust_eqtl_gene_count <-all_gene_count$eqtl_gene_count/all_gene_count$length*1000
all_gene_count1 <-filter(all_gene_count,cluster %in%c(1:7))

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

p1 <- ggplot(all_gene_count1, aes(x=as.factor(cluster), y=log(adjust_eqtl_gene_count),group=cluster,fill = as.factor(cluster))) + 
    geom_boxplot()+
    labs(x="Cluster",y="Log(number of eQTL-eGene pairs per kb)")+
    scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    theme_bw()+
    ggtitle("eQTL-eGene pairs")+
    theme(legend.position="none")+
    p_theme
pdf("13_eqtl_gene_count_boxplot.pdf",height=4.5,width=4.5)
print(p1)
dev.off()


p <- ggplot(all_gene_count1, aes(x=log(adjust_eqtl_gene_count),color=as.factor(cluster))) + 
geom_density(size=0.8)+
scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(number of eQTL-eGene pairs per kb)",color="Cluster") +
# coord_cartesian(xlim = c(0, 10))
ggtitle("eQTL-eGene pairs")+
p_theme
pdf("13_densityplot_eqtl_gene_count.pdf",height=4.5,width=4.9)
print(p)
dev.off()

