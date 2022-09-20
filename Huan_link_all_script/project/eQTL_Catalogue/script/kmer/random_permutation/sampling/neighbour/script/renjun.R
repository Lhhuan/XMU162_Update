load("09_kpca_rbfdot_chr1.Rdata")

pcn = 20
k = 50

# reduction
dat = pca1@pcv[,1:pcn]
umap = umap::umap(dat, method = "umap-learn")
tsne = Rtsne::Rtsne(dat, initial_dims = pcn, pca = FALSE)

# clustering
dat_knn = RcppHNSW::hnsw_knn(X = dat, k = k)
dat_edge = reshape2::melt(dat_knn$idx)
dat_graph = graph_from_edgelist(as.matrix(dat_edge[,c(1,3)]), directed = FALSE)
c = cluster_louvain(dat_graph)

# plot
df = data.frame(dat[,1:2], umap$layout, tsne$Y, as.factor(c$membership) )
colnames(df) = c(paste0("PC_",1:2), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
p <- ggplot(df, aes(x = PC_1, y = PC_2, color = cluster)) + 
  geom_point(size = 2.3, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave("pca.png", p, width = 7, height = 6)

p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 2.3, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave("tsne.png", p, width = 7, height = 6)

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 2.3, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave("umap.png", p, width = 7, height = 6)