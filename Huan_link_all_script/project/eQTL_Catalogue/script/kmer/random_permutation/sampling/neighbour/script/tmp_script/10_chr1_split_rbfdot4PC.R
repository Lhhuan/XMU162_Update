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
library(kernlab)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr")
load("09_kpca_rbfdot_chr1.Rdata")
load("../08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
Sorg_all <-Sorg
aa <- as.data.frame(pca1@pcv)
colnames(aa)[1:2] <- c("PC1","PC2")
pca <-aa

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 
#----------------------------
load("10_kpca_rbfdot_chr1_tsne.Rdata")
load("10_kpca_rbfdot_chr1_umap.Rdata")

#------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------kmeans
library(amap)
dat1 <-as.matrix(pca[,1:30])
kr_n <-function(i){
    print(paste(i,"start",sep="\t"))
    kr <-Kmeans(dat1, i)
    print(paste(i,"end",sep="\t"))
    return(kr)
}

re <-lapply(c(2:9),kr_n)

save(re,file="10_chr1_kmeans_2_9_30pca.Rdata")
plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}


plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_tsne_30PCA.png",p,dpi=300,width=10.7,height=4.3)
# #_-------------------------------------------------------------------------

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_UMAP_30PCA.png",p,dpi=150,width=10.7,height=4.3)

plot_pca <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(PC1=aa$PC1,PC2=aa$PC2,cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = PC1, y = PC2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}
plist <- lapply(c(1:8),plot_pca)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_Kpca_30PCA.png",p,dpi=150,width=10.7,height=4.3)





#-----------------------------------


#---------------------------------------------------------------------------------------------------------------------------------------kmeans
library(amap)
dat1 <-as.matrix(pca[,1:3])
kr_n <-function(i){
    print(paste(i,"start",sep="\t"))
    kr <-Kmeans(dat1, i)
    print(paste(i,"end",sep="\t"))
    return(kr)
}

re <-lapply(c(2:9),kr_n)

save(re,file="10_chr1_kmeans_2_9_3pca.Rdata")
plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}


plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_tsne_3PCA.png",p,dpi=300,width=10.7,height=4.3)
# #_-------------------------------------------------------------------------

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_UMAP_3PCA.png",p,dpi=150,width=10.7,height=4.3)

plot_pca <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(PC1=aa$PC1,PC2=aa$PC2,cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = PC1, y = PC2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}
plist <- lapply(c(1:8),plot_pca)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_Kpca_3PCA.png",p,dpi=150,width=10.7,height=4.3)

#----------------------------------------------------------------------------------------------------------------------------------------------------------


library(amap)
dat1 <-as.matrix(pca[,1:4])
kr_n <-function(i){
    print(paste(i,"start",sep="\t"))
    kr <-Kmeans(dat1, i)
    print(paste(i,"end",sep="\t"))
    return(kr)
}

re <-lapply(c(2:9),kr_n)

save(re,file="10_chr1_kmeans_2_9_4pca.Rdata")
plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_tsne_4PCA.png",p,dpi=300,width=10.7,height=4.3)
# #_-------------------------------------------------------------------------
plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_UMAP_4PCA.png",p,dpi=150,width=10.7,height=4.3)

plot_pca <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(PC1=aa$PC1,PC2=aa$PC2,cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = PC1, y = PC2, colour = factor(cluster)), size = 0.1, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}
plist <- lapply(c(1:8),plot_pca)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("10_chr1_kmeans_Kpca_4PCA.png",p,dpi=150,width=10.7,height=4.3)
#-------------------------------------------------------------------leiden
library(parallelDist)
sample_dat <-pca[,1:2]
eu = parallelDist::parDist(x = as.matrix(sample_dat), method = "euclidean",threads = 20)
cutoff <-quantile(eu,0.75) #22.17723
names(cutoff)=NULL

tu <-eu
tu[which(tu<=cutoff)]=0
tu[which(tu>cutoff)]=1
library(igraph)
# g1 <- graph_from_adjacency_matrix(as.matrix(tu),mode="undirected")
g1 <- graph_from_adjacency_matrix(as.matrix(eu),weight=TRUE,mode="undirected")
# g1 <- graph_from_adjacency_matrix(as.matrix(eu),mode="undirected")
# g <-cluster_leiden(g1)
# leiden_r <-data.frame(hotspot=rownames(sample_dat),cluster=g$membership)

# df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)

# p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                                                 panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
#                                                 axis.title.x = element_text(size = 10),
#                                                 axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

# p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
#     guides(col=guide_legend(title='Cluster'))
# ggsave("10_2PCA_chr1_euc_leiden_Umap_cluster.png",p,dpi=300,width=5,height=5)

#--------------------------------------------------cluster_louvain
g <-cluster_louvain(g1, weights = NULL)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("10_2PCA_chr1_euc_louvain_Umap_cluster.png",p,dpi=300,width=5,height=5)

df <-data.frame(tsne1=tsne$Y[,1],tsne2=tsne$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("10_2PCA_chr1_euc_louvain_tsne_cluster.png",p,dpi=300,width=5,height=4.5)


#--------------------------------------------------
#-------------------------------------------------------------------leiden
library(parallelDist)
sample_dat <-pca[,1:2]
eu = parallelDist::parDist(x = as.matrix(sample_dat), method = "euclidean",threads = 20)
cutoff <-quantile(eu,0.75) #22.17723
names(cutoff)=NULL

tu <-eu
tu[which(tu<=cutoff)]=0
tu[which(tu>cutoff)]=1
library(igraph)
# g1 <- graph_from_adjacency_matrix(as.matrix(tu),mode="undirected")
g1 <- graph_from_adjacency_matrix(as.matrix(eu),weight=TRUE,mode="undirected")
# g1 <- graph_from_adjacency_matrix(as.matrix(eu),mode="undirected")
# g <-cluster_leiden(g1)
# leiden_r <-data.frame(hotspot=rownames(sample_dat),cluster=g$membership)

# df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)

# p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                                                 panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
#                                                 axis.title.x = element_text(size = 10),
#                                                 axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

# p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
#     guides(col=guide_legend(title='Cluster'))
# ggsave("10_2PCA_chr1_euc_leiden_Umap_cluster.png",p,dpi=300,width=5,height=5)

#--------------------------------------------------cluster_louvain
g <-cluster_louvain(g1, weights = NULL)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("10_2PCA_chr1_euc_louvain_Umap_cluster.png",p,dpi=300,width=5,height=5)

df <-data.frame(tsne1=tsne$Y[,1],tsne2=tsne$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("10_2PCA_chr1_euc_louvain_tsne_cluster.png",p,dpi=300,width=5,height=4.5)


#--------------------------------------------------

















# #----------------------------------------------------------mcluster
# mc <-function(i){
#     print(paste(i,"start",sep="\t"))
#     mod <-Mclust(sample_dat,G=1:i)
#     print(paste(i,"end",sep="\t"))
#     return(mod)
# }
# models <-lapply(c(2:9),mc)
# plot <-function(i){
#     re <-data.frame(models[[i]]$classification)
#     colnames(re)[1]="cluster"
#     df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re$cluster)
#     p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.2, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
#     guides(col=guide_legend(title='Cluster'))
#     return(p)
# }
# plist <- lapply(c(1:8),plot)
# p =CombinePlots(plist,ncol=4,nrow=2)
# ggsave("10_2PCA_chr1_mclust_tsne.png",p,dpi=120,width=12,height=5)

# uplot <-function(i=NULL){
#     re <-data.frame(models[[i]]$classification)
#     colnames(re)[1]="cluster"
#     df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re$cluster)
#     p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.2, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
#     guides(col=guide_legend(title='Cluster'))
# }

# plist <- lapply(c(1:8),uplot)
# p =CombinePlots(plist,ncol=4,nrow=2)
# ggsave("10_2PCA_chr1_mclust_UMAP.png",p,dpi=120,width=12,height=5)

# #-------------------------------------------------
# library("kernlab")

# specc_n <-function(i){
#     print(paste(i,"start",sep="\t"))
#     sc <-specc(as.matrix(sample_dat),centers=i)
#     print(paste(i,"end",sep="\t"))
#     return(sc)
# }
# re <-lapply(c(2:9),specc_n)


# plot <-function(i){
#     print(paste(i,"start",sep="\t"))
#     df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]@.Data)
#     p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
#     guides(col=guide_legend(title='Cluster'))
#     print(paste(i,"end",sep="\t"))
#     return(p)
# }

# plist <- lapply(c(1:8),plot)
# p =CombinePlots(plist,ncol=4,nrow=2)
# ggsave("10_2PCA_chr1_specc_UMAP.png",p,dpi=150,width=12,height=5)

# plot <-function(i){
#     print(paste(i,"start",sep="\t"))
#     df <-data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2],cluster = re[[i]]@.Data)
#     p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
#     guides(col=guide_legend(title='Cluster'))
#     print(paste(i,"end",sep="\t"))
#     return(p)
# }
# plist <- lapply(c(1:8),plot)
# p =CombinePlots(plist,ncol=4,nrow=2)
# ggsave("10_2PCA_chr1_specc_tsne.png",p,dpi=150,width=12,height=5)