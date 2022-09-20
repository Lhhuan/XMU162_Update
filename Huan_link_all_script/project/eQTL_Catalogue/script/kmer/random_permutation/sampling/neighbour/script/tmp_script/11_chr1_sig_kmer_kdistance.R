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
library(mclust)
library(kmer)


setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr")
load("09_kpca_rbfdot_chr1.Rdata")
load("../08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
Sorg_all <-Sorg
aa <- as.data.frame(pca1@rotated)
colnames(aa)[1:2] <- c("PC1","PC2")
pca <-aa

chr="chr1:"
Sorg_chr1 <-Sorg_all[grep(chr,rownames(Sorg_all)),]


.kdist <- function(x, from, to, seqlengths, k) {
    .Call('_kmer_kdist', PACKAGE = 'kmer', x, from, to, seqlengths, k)
}


seqalongx=seq_along(1:nrow(Sorg_chr1))
k=6
kcounts <-as.matrix(Sorg_chr1)
seqlengths =apply(kcounts, 1, sum) + k - 1


d <- .kdist(kcounts, from = seqalongx - 1, to = seqalongx - 1,
                        seqlengths = seqlengths, k = k)

save(d,file ="11_kdistance_chr1.Rdata")
cutoff <-quantile(d,0.75,na.rm=TRUE)

names(cutoff)=NULL
tu <-d 
tu[which(is.na(tu))]=0
tu[which(tu<=cutoff)]=0
tu[which(tu>cutoff)]=1
library(igraph)
#-----------------------------------
g1 <- graph_from_adjacency_matrix(tu,mode="undirected")

load("10_kpca_rbfdot_chr1_tsne.Rdata")
load("10_kpca_rbfdot_chr1_umap.Rdata")

g <-cluster_leiden(g1, resolution_parameter=0.03,weights = NULL)
g_leiden<-g
g$membership[which(g$membership>=2)]=2
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("11_chr1_kdistance_leiden_Umap.png",p,dpi=300,width=5,height=4.5)

df <-data.frame(tsne1=tsne$Y[,1],tsne2=tsne$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("11_chr1_kdistance_leiden_tsne.png",p,dpi=300,width=5,height=4.5)


#----------------------------------------------------------------------------------cluster_louvain
g <-cluster_louvain(g1, weights = NULL)
save(g,file="11_chr1_kdistance_louvain.Rdata")
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("11_chr1_kdistance_kdistance_louvain_Umap.png",p,dpi=300,width=5,height=4.5)

df <-data.frame(tsne1=tsne$Y[,1],tsne2=tsne$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("11_chr1_kdistance_louvain_tsne.png",p,dpi=300,width=5,height=4.5)


df <-data.frame(PC1=aa$PC1,PC2=aa$PC2,cluster = g$membership)
p =ggplot(df) + geom_point(aes(x = PC1, y = PC2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("11_chr1_kdistance_louvain_kpca.png",p,dpi=300,width=5,height=4.5)

#--------------------------------

g1 <- graph_from_adjacency_matrix(as.matrix(d),weight=TRUE,,mode="undirected")

g <-cluster_louvain(g1, weights = NULL)
save(g,file="11_kmer_chr1_kdistance_louvain_all.Rdata")
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("11_kmer_chr1_kdistance_kdistance_louvain_Umap_all.png",p,dpi=300,width=5,height=4.5)

df <-data.frame(tsne1=tsne$Y[,1],tsne2=tsne$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("11_kmer_chr1_kdistance_louvain_tsne_all.png",p,dpi=300,width=5,height=4.5)


df <-data.frame(PC1=aa$PC1,PC2=aa$PC2,cluster = g$membership)
p =ggplot(df) + geom_point(aes(x = PC1, y = PC2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("11_kmer_chr1_kdistance_louvain_kpca_all.png",p,dpi=300,width=5,height=4.5)