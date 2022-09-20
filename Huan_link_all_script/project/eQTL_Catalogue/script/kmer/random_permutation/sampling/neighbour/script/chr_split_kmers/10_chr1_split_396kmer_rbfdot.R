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
pcn = 30


# reduction
dat = pca1@pcv[,1:pcn]
# umap = umap::umap(dat, method = "umap-learn")
set.seed(123)
umap <- uwot::umap(dat)
set.seed(123)
tsne = Rtsne::Rtsne(dat, initial_dims = pcn, pca = FALSE)

load("../08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
Sorg_all <-Sorg


chr="chr1:"
Sorg_chr1 <-Sorg_all[grep(chr,rownames(Sorg_all)),]
# clustering
k = 5
dat_knn = RcppHNSW::hnsw_knn(X = as.matrix(Sorg_chr1), k = k)
dat_edge = reshape2::melt(dat_knn$idx)
dat_graph = igraph::graph_from_edgelist(as.matrix(dat_edge[,c(1,3)]), directed = FALSE)
c = igraph::cluster_louvain(dat_graph)

# plot
# df = data.frame(dat[,1:2], umap$layout, tsne$Y, as.factor(c$membership) )
df = data.frame(dat[,1:2], umap, tsne$Y, as.factor(c$membership) )
colnames(df) = c(paste0("PC_",1:2), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
df$hotspot <- rownames(as.data.frame(pca1@xmatrix))

p <- ggplot(df, aes(x = PC_1, y = PC_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_365kmer_louvain_pca",pcn,"_k",k,"_pca.png"), p, width = 6, height = 5.3)

p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_365kmer_louvain_pca",pcn,"_k",k,"_tsne.png"), p, width = 6, height = 5.3)

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 2, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("10_chr1_365kmer_louvain_pca",pcn,"_k",k,"_umap.png"), p, width = 6, height = 5.3)

# write.table(df,paste0("10_chr1_365kmer_louvain_pca",pcn,"_k",k,".txt"),row.names = F, col.names = T,quote =F,sep="\t")