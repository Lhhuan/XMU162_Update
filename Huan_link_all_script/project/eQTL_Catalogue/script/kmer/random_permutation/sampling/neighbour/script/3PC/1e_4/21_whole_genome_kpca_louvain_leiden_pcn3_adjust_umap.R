setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/")
load("whole_genome_rbfdot_pca.Rdata")
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/UMAP_adjust/")

library(ggplot2)
library(igraph)
library(dplyr)
library(ggsci)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                #   axis.title.y = element_text(size = 10),
                                                #   axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))

pcn = 3

# reduction
dat1 = pca1@pcv[,1:pcn]
rownames(dat1)=1:nrow(dat1)
dat_dup <-dat1[duplicated(dat1),]
dat <-dat1[!duplicated(dat1),]

k = 50
set.seed(1)
dat_knn = RcppHNSW::hnsw_knn(X = dat, k = k)
dat_edge = reshape2::melt(dat_knn$idx)
dat_graph = igraph::graph_from_edgelist(as.matrix(dat_edge[,c(1,3)]), directed = FALSE)
resolution=0.0001
set.seed(40)
g <-cluster_leiden(dat_graph, resolution_parameter=resolution)

dist = seq(0.005,0.2,0.001)
# neighbor = seq(10,20,1)
neighbor = seq(12,20,1)
for(i in neighbor){
  for(j in dist){
  set.seed(123)
  umap <- uwot::umap(dat,n_neighbors = i,min_dist = j)
  df = data.frame(dat[,1:3], umap, as.factor(g$membership) )
  colnames(df) = c(paste0("PC_",1:3), paste0("UMAP_",1:2), "cluster")
  df$hotspot <- rownames(as.data.frame(pca1@xmatrix))[as.numeric(rownames(dat))]
  #-----------------------
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
  ggsave(paste0("21_neig_",i,"_mindist_",j,"_umap_leiden_pca3_k50_resolution_1e_4.png"), p, width = 5, height = 4.3)
  print(c(i,j))
  }
}

# # c = igraph::cluster_louvain(dat_graph)



# df = data.frame(dat[,1:3], umap, tsne$Y, as.factor(g$membership) )
# # df = data.frame(dat[,1:4], umap, tsne$Y, as.factor(g$membership) )
# colnames(df) = c(paste0("PC_",1:3), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
# df$hotspot <- rownames(as.data.frame(pca1@xmatrix))[as.numeric(rownames(dat))]


# p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
#   geom_point(size = 0.8, alpha = 0.6) + 
#   scale_color_d3("category20") + 
#   theme_bw()
# ggsave(paste0("11_whole_genome_umap_leiden_pca",pcn,"_k",k,"_resolution",resolution,"_tsne.png"), p, width = 6, height = 5.3)

# p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
#   geom_point(size = 0.8, alpha = 0.6) + 
#   scale_color_d3("category20") + 
#   theme_bw()
# ggsave(paste0("11_whole_genome_umap_leiden_pca",pcn,"_k",k,"_resolution",resolution,"_umap.png"), p, width = 6, height = 5.3)

# dup <-data.frame(dat_dup[,1:3],hotspot=rownames(as.data.frame(pca1@xmatrix))[as.numeric(rownames(dat_dup))])
# colnames(dup)[1:3] <-paste0("PC_",1:3)
# dup <-left_join(dup,df[,1:8],by=c("PC_1","PC_2","PC_3"))
# dup <-dup[,c(1:3,5:9,4)]
# fdatl <- bind_rows(dup,df)

# write.table(fdatl,paste0("11_whole_genome_umap_leiden_pca",pcn,"_k",k,"_resolution",resolution,".txt"),row.names = F, col.names = T,quote =F,sep="\t")
