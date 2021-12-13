library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(Seurat)
library(R.utils)
library(reshape2)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

#---------------------------------------------------------
library(Hmisc)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")

org<-read.table("../../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$chr,org$start,org$end,sep="_")
org <-org[,c("Marker","overlap_fraction","hotspot")]
org1 <-org%>%group_by(Marker,hotspot)%>%summarise(sum_overlap_fraction =sum(overlap_fraction))%>%data.frame()
org2 <-dcast(org1[,c(2,1,3)],hotspot~Marker)

GC<-read.table("../../../output/Tissue_merge/Cis_eQTL/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833_GC.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(GC)[1:5] <-c("chr","start","end","AT_content","GC_content")
GC$hotspot <-paste(GC$chr,GC$start,GC$end,sep="_")
GC <-GC[,c("hotspot","GC_content")]

dat <- inner_join(GC,org2,by="hotspot")
rownames(dat) = dat[,"hotspot"]
dat = dat[,-1]

save(dat,file="04_hotspot_mark_fraction.Rdata")
library(irlba)
set.seed(11)
pca = prcomp_irlba(dat, n = 5)
pca_result = pca$x
rownames(pca_result) = rownames(dat)
df <-data.frame(pca_result) 
p =ggplot(df) + geom_point(aes(x = PC1, y = PC2), size = 0.8, alpha = 0.3)+p_theme #, colour = "#6a737b"
ggsave("04_pca_cluster_GC.png",p,dpi=300,width=7,height=7)

#---------------

dat=dat[,-1] #-GC_content

set.seed(11)
pca = prcomp_irlba(dat, n =4)
pca_result = pca$x
rownames(pca_result) = rownames(dat)
df <-data.frame(pca_result) 
p =ggplot(df) + geom_point(aes(x = PC1, y = PC2), size = 0.8, alpha = 0.3)+p_theme #, colour = "#6a737b"
ggsave("04_pca_cluster.png",p,dpi=300,width=7,height=7)

#------------------------------------

eu_dist <- dist(dat, method = "euclidean", diag = T, upper = T, p = 2)
save(eu_dist,file="04_10_markers_euclidean_dist.Rdata")