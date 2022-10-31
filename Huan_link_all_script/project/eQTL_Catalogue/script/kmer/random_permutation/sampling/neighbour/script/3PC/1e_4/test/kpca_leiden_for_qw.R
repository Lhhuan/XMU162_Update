library(unix)
rlimit_as(1e18)
library(kernlab)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure")
load("08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")

#kpca
pca1<- kpca(~.,data=Sorg, kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)
save(pca1,file="whole_genome_rbfdot_pca.Rdata")

#=========================
pcn = 3
# reduction
dat = pca1@rotated[,1:pcn]
rownames(dat)=1:nrow(dat)

#tsne and umap
set.seed(123)
umap <- uwot::umap(dat)
save(umap,file=paste0("11_whole_genome_",pcn,"kpca_umap.Rdata"))
set.seed(123)
tsne = Rtsne::Rtsne(dat, initial_dims = pcn, pca = FALSE)
save(tsne,file=paste0("11_whole_genome_",pcn,"kpca_tsne.Rdata"))

# graph
k = 50
set.seed(1)
dat_knn = RcppHNSW::hnsw_knn(X = dat, k = k)
dat_edge = reshape2::melt(dat_knn$idx)
dat_graph = igraph::graph_from_edgelist(as.matrix(dat_edge[,c(1,3)]), directed = FALSE)
#clustering       leiden
resolution=0.0001
set.seed(40)#k=30
g <-cluster_leiden(dat_graph, resolution_parameter=resolution)
#plot
df = data.frame(dat, umap, tsne$Y, as.factor(g$membership) )
colnames(df) = c(paste0("PC_",1:3), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
#----
p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("11_whole_genome_leiden_pca",pcn,"_k",k,"_resolution",resolution,"_tsne.png"), p, width = 6, height = 5.3)
#-----
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("11_whole_genome_leiden_pca",pcn,"_k",k,"_resolution",resolution,"_umap.png"), p, width = 6, height = 5.3)