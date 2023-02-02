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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/Cis_Regulatory_Elements/")
org <-read.table("hotspot_cluster_cCREs.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)<-c("CHR","start","end","cluster","cre_chr","cre_start","cre_end","name","score","strand","thickStart","thickEnd","reserved","ccre","encodeLabel","zScore","ucscLabel","accessionLabel","description")
org$hotspot <-paste(org$CHR,org$start,org$end,sep=":")
org$cluster <-as.factor(org$cluster)
org <- org[,c("hotspot","cluster","ucscLabel","name")]
org1 <- unique(org[,c("hotspot","cluster","ucscLabel")])

dat <-org1%>%group_by(ucscLabel,cluster)%>%summarise(ucscLabel_count=n())%>%data.frame()

all_hotspot <- read.table("../GWAS/17_whole_genome_leiden_pca3_k50_resolution_2e-5_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(all_hotspot) <-c("CHR","start","end","cluster")
all_hotspot$cluster <-as.factor(all_hotspot$cluster)
all_hotspot$hotspot <- paste(all_hotspot$CHR,all_hotspot$start,all_hotspot$end,sep=":")
all_hotspot <- all_hotspot[,c("hotspot","cluster")]
cluster_n <- all_hotspot %>%group_by(cluster)%>%summarise(hotspot_count=n())%>%data.frame()


dat <- left_join(dat,cluster_n,by=c("cluster"))
dat$proporation <-dat$ucscLabel_count/dat$hotspot_count

fdat <-dcast(dat[,c(1,2,5),],cluster~ucscLabel)

color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
# color2 = colorRampPalette(c('#D3BDF3','#9469FE','#5C33FC'))(50)
library(pheatmap)
pdf("22_2_cCRE_propotion_heatmap.pdf",width=2.7,height=1.2)
pheatmap(fdat[,2:6],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = TRUE,show_rownames = T,show_colnames = T, cellwidth =27, cellheight = 27)
dev.off()


