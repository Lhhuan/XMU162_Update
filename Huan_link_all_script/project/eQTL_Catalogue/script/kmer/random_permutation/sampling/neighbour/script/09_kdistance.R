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

.kdist <- function(x, from, to, seqlengths, k) {
    .Call('_kmer_kdist', PACKAGE = 'kmer', x, from, to, seqlengths, k)
}

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/")
load("08_random_5000_in_permutation_wilocx_overlap_sig.Rdata")
Sorg <-random5000
seqalongx=seq_along(1:nrow(Sorg))
k=6
kcounts <-as.matrix(Sorg)
seqlengths =apply(kcounts, 1, sum) + k - 1


d <- .kdist(kcounts, from = seqalongx - 1, to = seqalongx - 1,
                        seqlengths = seqlengths, k = k)

save(d,file ="09_kdistance_random5000.Rdata")

cutoff <-quantile(d,0.75)

names(cutoff)=NULL
tu <-d 
tu[which(tu<=cutoff)]=0
tu[which(tu>cutoff)]=1
library(igraph)
#-----------------------------------
g1 <- graph_from_adjacency_matrix(tu,mode="undirected")
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/")
load("08_new_pca_tsne.Rdata")
load("08_new_pca_umap.Rdata")

g <-cluster_leiden(g1, resolution_parameter=0.03,weights = NULL)
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("09_random5000_leiden_Umap_cluster_new.png",p,dpi=300,width=5,height=4.5)

df <-data.frame(tsne1=tsne$Y[,1],tsne2=tsne$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("09_random5000_leiden_tsne_cluster_new.png",p,dpi=300,width=5,height=4.5)


#----------------------------------------------------------------------------------cluster_louvain
g <-cluster_louvain(g1, weights = NULL)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("09_random5000_kdistance_louvain_Umap_cluster_new.png",p,dpi=300,width=5,height=5)

df <-data.frame(tsne1=tsne$Y[,1],tsne2=tsne$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("09_random5000_louvain_tsne_cluster_new.png",p,dpi=300,width=5,height=4.5)