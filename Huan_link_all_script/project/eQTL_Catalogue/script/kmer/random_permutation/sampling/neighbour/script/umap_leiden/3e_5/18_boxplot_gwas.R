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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/re_3e_5/GWAS/")
org <-read.table("17_whole_genome_umap_leiden_pca5_k50_resolution3e-05_gwas_metaid.bed.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$CHR,org$start,org$end,sep=":")
org$length <- org$end-org$start

org$cluster <-as.factor(org$cluster)
org$trait_gwas <- paste(org$GWAS_chr,org$GWAS_start,org$GWAS_end,org$meta_id,sep="_")
dat <-org%>%group_by(hotspot,length,cluster)%>%summarise(trait_count=n())%>%data.frame()
dat$adjust_trait_count <-dat$trait_count/dat$length*1000

dat$cluster <-as.factor(dat$cluster)

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))


p1 <- ggplot(dat, aes(x=as.factor(cluster), y=log(adjust_trait_count),group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
theme_bw()+
labs(x="Cluster",y="Log(number of snp-trait pairs per kb)")+
ggtitle("SNP-trait pairs")+
theme(legend.position ="none")+
p_theme

pdf("18_boxpolt_whole_genome_gwas_trait.pdf",height=4.7,width=4.5)
print(p1)
dev.off()



p <- ggplot(dat, aes(x=log(adjust_trait_count),color=cluster)) + 
geom_density(size=0.8)+
scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(Number of snp-trait pairs per kb)",color="Cluster") +
ggtitle("SNP-trait pairs")+
p_theme

pdf("18_densityplot_whole_genome_gwas_trait.pdf",height=4.5,width=4.9)
print(p)
dev.off()

