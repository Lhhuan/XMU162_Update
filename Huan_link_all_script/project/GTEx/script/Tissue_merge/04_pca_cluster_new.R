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
library(umap)
library(uwot)
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

#---------------------------------------------------------
library(Hmisc)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")

org<-read.table("../../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$chr,org$start,org$end,sep="_")
org <-org[,c("Marker","overlap_fraction","hotspot")]
org1 <-org%>%group_by(Marker,hotspot)%>%summarise(sum_overlap_fraction =sum(overlap_fraction))%>%data.frame()
org2 <-dcast(org1[,c(2,1,3)],hotspot~Marker)
rownames(org2)<-org2$hotspot
dat <-org2[,-1]
dat1 <-t(dat)
pca=prcomp(dat1)
# dat1 <-dat1[-1,]
# pca_s <-RunPCA(dat1,npcs=50)

# library(irlba)
# set.seed(11)
# pca = prcomp_irlba(dat, n = 9)
pca_result = pca$rotation
df <-data.frame(pca_result) 
p =ggplot(df) + geom_point(aes(x = PC1, y = PC2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("04_pca_cluster_new.png",p,dpi=150,width=5,height=5)

tsne = Rtsne::Rtsne(pca_result, pca = FALSE,check_duplicates = FALSE)
save(tsne,file="04_new_pca_tsne.Rdata")

df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("04_PCA_tsn_cluster_new.png",p,dpi=150,width=5,height=5)
#---------------
rumap <- uwot::umap(pca_result)
rownames(rumap)=rownames(dat)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2])
save(rumap,file="04_new_pca_umap.Rdata")

p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("04_PCA_Umap_cluster_new.png",p,dpi=150,width=5,height=5)
#------------------------------------

# eu_dist <- dist(dat, method = "euclidean",diag = TRUE,upper = T)
# # g <-graph_from_adjacency_matrix(eu_dist)
# # df <- melt(as.matrix(eu_dist), varnames = c("hotspot1", "hotspot2"))
# # colnames(df)[3] <-"euclidean_dist"
# # save(df,file="04_10_markers_euclidean_dist_df.Rdata")
# save(eu_dist,file="04_10_markers_euclidean_dist_matrix.Rdata")



#--------------------


# pca = prcomp_irlba(dat,n=50)
# pca1 <-data.frame(pca)
# impPrComp = princomp(dat)
# screeplot(impPrComp,type="line")