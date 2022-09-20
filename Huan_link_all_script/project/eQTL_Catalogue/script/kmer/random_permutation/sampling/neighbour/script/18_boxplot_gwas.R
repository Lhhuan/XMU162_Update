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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/GWAS/")
org <-read.table("17_chr1_louvain_pca5_k500_hotspot_gwas_metaid.bed.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$CHR,org$start,org$end,sep=":")

org$cluster <-as.factor(org$cluster)
dat <-org%>%group_by(hotspot,cluster)%>%summarise(trait_count=n())%>%data.frame()
dat$cluster <-as.factor(dat$cluster)
p1 <- ggplot(dat, aes(x=as.factor(cluster), y=trait_count,group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
theme_bw()+
theme(legend.position="none")+
labs(x="Cluster",y="Number of snp-trait pairs")+
pdf("18_boxpolt_chr1_gwas_trait.pdf",height=6,width=6)
print(p1)
dev.off()



p <- ggplot(dat, aes(x=trait_count,color=cluster)) + 
geom_density()+
scale_color_d3("category20") + 
theme_bw()+
xlab("Number of snp-trait pairs")+
coord_cartesian(xlim = c(0, 50))

pdf("18_densityplot_chr1_gwas_trait.pdf",height=6,width=7)
print(p)
dev.off()



