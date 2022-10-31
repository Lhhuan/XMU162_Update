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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/3pca_5e_6/CPG_island/")
org <- read.table("../../../11_whole_genome_umap_leiden_pca3_k30_resolution5e-06.txt",header = T,sep = "\t") %>% as.data.frame()
# org <-read.table("17_whole_genome_leiden_pca5_k30_resolution1e-04_gwas_metaid.bed.gz",header = T,sep = "\t") %>% as.data.frame()
cpg <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/CPG_island_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()


org$cluster <-as.factor(org$cluster)

colnames(cpg)[1:3] <-c("CHR","start","end")
# cpg$hotspot <-paste(cpg$CHR,cpg$start,cpg$end,sep=":")
cpg$hotspot <-paste0(cpg$CHR,":",cpg$start,"-",cpg$end)
cpg <-cpg[,c("CHR","hotspot")]
cpg <-left_join(cpg,org[,c("hotspot","cluster")],by="hotspot")



all_dat <-org%>%group_by(cluster)%>%summarise(cluster_count=n())%>%data.frame()
cpg_dat <- cpg%>%group_by(cluster)%>%summarise(cpg_cluster_count=n())%>%data.frame()
dat <-left_join(cpg_dat,all_dat,by="cluster")
dat$proportion <- dat$cpg_cluster_count/dat$cluster_count
dat1 <- filter(dat,cluster %in%c(1:7))

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))


p1 <- ggplot(dat1, aes(x=as.factor(cluster), y=proportion,group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_bar(stat = 'identity', width=0.8)+
scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1"))+
theme_bw()+
labs(x="Cluster",y="Fraction of hotspots")+
ggtitle("CPG island")+
theme(legend.position ="none")+
p_theme

pdf("19_barpolt_whole_genome_cpg_island.pdf",height=4.7,width=4.5)
print(p1)
dev.off()



