setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr")
load("09_kpca_rbfdot_chr1.Rdata")
library(ggplot2)
library(igraph)
library(ggsci)
variances <- data.frame(pca1@eig)
variances$Dimensions <- 1:nrow(variances)
variances$y <- variances$pca1.eig/sum(variances$pca1.eig) *100
variances$sum_y <-cumsum(variances$pca1.eig)/sum(variances$pca1.eig)*100
head(which(variances$sum_y>60))
pcn = 10


# reduction
dat = pca1@pcv[,1:pcn]
# umap = umap::umap(dat, method = "umap-learn")
set.seed(123)
umap <- uwot::umap(dat)
set.seed(123)
tsne = Rtsne::Rtsne(dat, initial_dims = pcn, pca = FALSE)

# clustering
k = 300
# dat_knn = RcppHNSW::hnsw_knn(X = dat, k = k)
dat_knn = RcppHNSW::hnsw_knn(X = umap, k = k)
dat_edge = reshape2::melt(dat_knn$idx)
dat_graph = igraph::graph_from_edgelist(as.matrix(dat_edge[,c(1,3)]), directed = FALSE)
c = igraph::cluster_louvain(dat_graph)

# plot
# df = data.frame(dat[,1:2], umap$layout, tsne$Y, as.factor(c$membership) )
df = data.frame(dat[,1:2], umap, tsne$Y, as.factor(c$membership) )
colnames(df) = c(paste0("PC_",1:2), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
p <- ggplot(df, aes(x = PC_1, y = PC_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_louvain_renjun_pca",pcn,"_k",k,"_pca.png"), p, width = 6, height = 5.3)

p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_louvain_renjun_pca",pcn,"_k",k,"_tsne.png"), p, width = 6, height = 5.3)

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_louvain_renjun_pca",pcn,"_k",k,"_umap.png"), p, width = 6, height = 5.3)

#----------
df = data.frame(dat[,1:2], umap, tsne$Y, as.factor(c$membership) )
colnames(df) = c(paste0("PC_",1:2), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")

p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_louvain_renjun_pca",pcn,"_k",k,"_cluster_umap_tsne.png"), p, width = 6, height = 5.3)

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_louvain_renjun_pca",pcn,"_k",k,"_cluster_umap_umap.png"), p, width = 6, height = 5.3)

#------------------

#leiden
resolution=0.019
g <-cluster_leiden(dat_graph, resolution_parameter=resolution)
df1 = data.frame(dat[,1:2], umap, tsne$Y, as.factor(g$membership) )
colnames(df1) = c(paste0("PC_",1:2), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
p <- ggplot(df1, aes(x = PC_1, y = PC_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_leiden_renjun_pca",pcn,"_k",k,"_resolution",resolution,"_pca.png"), p, width = 6, height = 5.3)

p <- ggplot(df1, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_leiden_renjun_pca",pcn,"_k",k,"_resolution",resolution,"_tsne.png"), p, width = 6, height = 5.3)

p <- ggplot(df1, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_leiden_renjun_pca",pcn,"_k",k,"_resolution",resolution,"_umap.png"), p, width = 6, height = 5.3)