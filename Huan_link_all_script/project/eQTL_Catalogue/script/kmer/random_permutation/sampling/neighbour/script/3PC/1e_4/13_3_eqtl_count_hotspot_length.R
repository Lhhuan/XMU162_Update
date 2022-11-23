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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/gene/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
Cluster <-read.table("../../11_whole_genome_leiden_pca3_k50_resolution1e-04.txt",header = T,sep = "\t") %>% as.data.frame()
tmp1 <-strsplit(Cluster$hotspot,":")
tmp2 <- do.call(rbind,tmp1)%>%data.frame()
Cluster$CHR <- tmp2$X1
cluster_chr <-Cluster %>%group_by(cluster,CHR)%>%summarise(count=n())%>%data.frame()
cluster_n <- Cluster%>%group_by(cluster)%>%summarise(cluster_count=n())%>%data.frame()
cluster_chr <-left_join(cluster_chr,cluster_n,by="cluster")
cluster_chr$ratio <- cluster_chr$count/cluster_chr$cluster_count


all_gene$hotspot <-paste0(all_gene$V1,":",all_gene$V2,"-",all_gene$V3)
all_gene$eqtl <- paste0(all_gene$V4,":",all_gene$V5,"-",all_gene$V6)
all_gene$length <- all_gene$V3 - all_gene$V2
colnames(all_gene)[7] <-"ENSG"
eqtl <- all_gene[,c("hotspot","eqtl","length")]%>%unique()
eqtl<-left_join(eqtl,Cluster[,c("cluster","hotspot")],by="hotspot")
eqtl_n <-eqtl%>%group_by(hotspot,cluster,length)%>%summarise(eqtl_count=n())%>%data.frame()
eqtl_n$adjust_eqtl <- eqtl_n$eqtl_count/eqtl_n$length*1000

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

eqtl_n$cluster <-factor(eqtl_n$cluster,levels=c(4,2,1,5,6,3))
p1 <- ggplot(eqtl_n, aes(x=cluster, y=log(adjust_eqtl),fill =cluster)) + 
    geom_boxplot()+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
    scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    theme_bw()+
    ggtitle("eQTL")+
    p_theme +
    # labs(x="Cluster",y="mean(signal of hic)")
    labs(x="Cluster",y="Log(Number of eQTL per kb)")
pdf("13_3_boxplot_of_adjust_eqtl_count_cluster.pdf",height=3.1,width=3)
print(p1)
dev.off() 

#====================
f4 <-filter(eqtl_n,cluster==4)
p <- ggplot(f4, aes(x=log(adjust_eqtl))) + 
geom_density(size=0.8,color="#C22324")+
theme_bw()+
labs(x="Log(Number of eQTL per kb)") +
ggtitle("eQTL(C4)")+
p_theme

pdf("13_3_densityplot_C4_eQTL.pdf",height=4.5,width=4.9)
print(p)
dev.off()

#===========================
eQTL_egene <- all_gene[,c("hotspot","eqtl","ENSG")] %>%unique()
eQTL_egene_n <- eQTL_egene%>%group_by(hotspot,eqtl)%>%summarise(ENSG_n =n())%>%data.frame()
eQTL_egene_n <- left_join(eQTL_egene_n,Cluster[,c("cluster","hotspot")],by="hotspot")

eQTL_egene_n$cluster <-factor(eQTL_egene_n$cluster,levels=c(4,2,1,5,6,3))

p1 <- ggplot(eQTL_egene_n, aes(x=cluster, y=log(ENSG_n),fill =cluster)) + 
    geom_boxplot()+
    # geom_violin()+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
    scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    theme_bw()+
    ggtitle("eQTL")+
    p_theme +
    # labs(x="Cluster",y="mean(signal of hic)")
    labs(x="Cluster",y="log(Number of egene per eQTL)")
pdf("13_3_boxplot_of_Number_of_egene_per_eQTL.pdf",height=4.7,width=4.5)
print(p1)
dev.off() 
#============================egene affected by no eqtl
egene_n <- eQTL_egene%>%group_by(hotspot,ENSG)%>%summarise(eqtl_n=n())%>%data.frame()
egene_n <- left_join(egene_n,Cluster[,c("cluster","hotspot")],by="hotspot")

egene_n$cluster <-factor(egene_n$cluster,levels=c(4,2,1,5,6,3))

p1 <- ggplot(egene_n, aes(x=cluster, y=log(eqtl_n),fill =cluster)) + 
    geom_boxplot()+
    # geom_violin()+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
    scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    theme_bw()+
    ggtitle("Number of eQTL per egene")+
    p_theme +
    # labs(x="Cluster",y="mean(signal of hic)")
    labs(x="Cluster",y="log(Number of eQTL per egene)")
pdf("13_3_boxplot_of_Number_of_eQTL_per_egene.pdf",height=4.7,width=4.5)
print(p1)
dev.off() 
#=================================NO. of egene in hotspot

eQTL_egene1 <- eQTL_egene[,c(1,3)]%>%unique()
egene_n1 <- eQTL_egene1%>%group_by(hotspot)%>%summarise(egene_n=n())%>%data.frame()
egene_n1 <- left_join(egene_n1,Cluster[,c("cluster","hotspot")],by="hotspot")
egene_n1 <-left_join(egene_n1,all_gene[,c("hotspot","length")],by="hotspot")
egene_n1$cluster <-factor(egene_n1$cluster,levels=c(4,2,1,5,6,3))
egene_n1$adjust_egene_n <- egene_n1$egene_n/egene_n1$length*1000
p1 <- ggplot(egene_n1, aes(x=cluster, y=log(adjust_egene_n),fill =cluster)) + 
    geom_boxplot()+
    # geom_violin()+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
    scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    theme_bw()+
    ggtitle("eGene")+
    p_theme +
    # labs(x="Cluster",y="mean(signal of hic)")
    labs(x="Cluster",y="Log(Number of egene per kb)")
pdf("13_3_boxplot_of_Number_of_egene_per_kb.pdf",height=3.1,width=3)
print(p1)
dev.off() 
#==================================================================

p1 <- ggplot(eqtl_n, aes(x=cluster, y=log10(length),fill =cluster)) + 
    geom_boxplot()+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_manual(values=c(1="#1E77B4",2="#FF7F0E",3="#2CA02C",4="#C22324",5="#9567BD",6="#8C554B"))+
    scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    theme_bw()+
    ggtitle("Length")+
    p_theme +
    # labs(x="Cluster",y="mean(signal of hic)")
    labs(x="Cluster",y="Log10(length of hotspot)")+
    scale_y_continuous(name="Length of hotspot(bp)", breaks=c(2,3,4,5), labels=c(100,1000,10000,100000))+
pdf("13_3_boxplot_of_length_of_hotspot_cluster.pdf",height=3.1,width=3.2)
print(p1)
dev.off() 