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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/GWAS/")
org <-read.table("17_whole_genome_leiden_pca3_k50_resolution1e-04_gwas_metaid.bed.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$CHR,org$start,org$end,sep=":")
org$length <- org$end-org$start

org$cluster <-as.factor(org$cluster)
org$trait_gwas <- paste(org$GWAS_chr,org$GWAS_start,org$GWAS_end,org$meta_id,sep="_")
dat <-org%>%group_by(hotspot,length,cluster)%>%summarise(trait_count=n())%>%data.frame()

all_hotspot <- read.table("17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(all_hotspot) <-c("CHR","start","end","cluster")
all_hotspot$cluster <-as.factor(all_hotspot$cluster)
all_hotspot$hotspot <- paste(all_hotspot$CHR,all_hotspot$start,all_hotspot$end,sep=":")
all_hotspot$length <- all_hotspot$end-all_hotspot$start
all_hotspot <- all_hotspot[,c("hotspot","length","cluster")]
all_hotspot <- left_join(all_hotspot,dat,by=c("hotspot","length","cluster"))
all_hotspot$trait_count[is.na(all_hotspot$trait_count)] <-0
all_hotspot$adjust_trait_count<-all_hotspot$trait_count/all_hotspot$length*1000

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

dat1 <-filter(all_hotspot,cluster %in%c(1:6))
dat1$cluster <-factor(dat1$cluster,levels=c(6,4,1,2,5,3))
p1 <- ggplot(dat1, aes(x=as.factor(cluster), y=log(adjust_trait_count),group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Log(number of snp-trait pairs per kb)")+
ggtitle("SNP-trait pairs")+
theme(legend.position ="none")+
p_theme

pdf("18_boxpolt_whole_genome_gwas_trait.pdf",height=4.7,width=4.5)
print(p1)
dev.off()



p1 <- ggplot(dat1, aes(x=as.factor(cluster), y=log(adjust_trait_count+0.1),group=cluster,fill = as.factor(cluster))) + 
# geom_boxplot(outlier.colour = NA)+
geom_boxplot()+
# scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
theme_bw()+
labs(x="Cluster",y="Log(number of snp-trait pairs per kb)")+
ggtitle("SNP-trait pairs")+
theme(legend.position ="none")+
p_theme

pdf("18_boxpolt_whole_genome_gwas_trait_1E-8.pdf",height=4.7,width=4.5)
print(p1)
dev.off()



p <- ggplot(dat1, aes(x=log(adjust_trait_count),color=cluster)) + 
geom_density(size=0.8)+
scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(Number of snp-trait pairs per kb)",color="Cluster") +
ggtitle("SNP-trait pairs")+
p_theme

pdf("18_densityplot_whole_genome_gwas_trait.pdf",height=4.5,width=4.9)
print(p)
dev.off()

all_dat <-dat1%>%group_by(cluster)%>%summarise(cluster_count=n())%>%data.frame()
aa <-filter(dat1,adjust_trait_count>0)
aa_dat <-aa%>%group_by(cluster)%>%summarise(GWAS_anno_count=n())%>%data.frame()
fdatcount <- left_join(all_dat,aa_dat,by="cluster")
fdatcount$propotion <- fdatcount$GWAS_anno_count /fdatcount$cluster_count