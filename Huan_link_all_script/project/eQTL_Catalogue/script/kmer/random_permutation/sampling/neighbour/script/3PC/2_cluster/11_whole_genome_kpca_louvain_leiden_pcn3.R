setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/")
load("whole_genome_rbfdot_pca.Rdata")
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome")

library(ggplot2)
library(igraph)
library(dplyr)
library(ggsci)
variances <- data.frame(pca1@eig)
variances$Dimensions <- 1:nrow(variances)
variances$y <- variances$pca1.eig/sum(variances$pca1.eig) *100
variances$sum_y <-cumsum(variances$pca1.eig)/sum(variances$pca1.eig)*100
head(which(variances$sum_y>60))

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                #   axis.title.y = element_text(size = 10),
                                                #   axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))


p1 <-ggplot(data = variances[1:30,], mapping = aes(x =Dimensions, y = y)) + 
    geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.9)+
    p_theme+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(breaks = seq(1, 30, by = 1))+
    labs(y="Percentage of explained variances",title="Scree plot")
ggsave(paste0("11_kPCA_screeplot_whole_genome.pdf"),p1,dpi=300,width=6,height=5) 

pcn = 3


# reduction
# dat1 = pca1@pcv[,1:pcn]
dat1 = pca1@rotated[,1:pcn]
# rownames(dat1)=1:nrow(dat1)
dat_dup <-dat1[duplicated(dat1),]
dat <-dat1[!duplicated(dat1),]

# umap = umap::umap(dat, method = "umap-learn")
set.seed(123)
umap <- uwot::umap(dat)
save(umap,file=paste0("11_whole_genome_",pcn,"kpca_umap.Rdata"))
set.seed(123)
tsne = Rtsne::Rtsne(dat, initial_dims = pcn, pca = FALSE)
save(tsne,file=paste0("11_whole_genome_",pcn,"kpca_tsne.Rdata"))
# clustering
k = 50
set.seed(1)
dat_knn = RcppHNSW::hnsw_knn(X = dat, k = k)
dat_edge = reshape2::melt(dat_knn$idx)
dat_graph = igraph::graph_from_edgelist(as.matrix(dat_edge[,c(1,3)]), directed = FALSE)
#==========================================
#                   leiden
resolution=0.00002
set.seed(40)#k=30
g <-cluster_leiden(dat_graph, resolution_parameter=resolution)

df = data.frame(dat[,1:3], umap, tsne$Y, as.factor(g$membership) )
# df = data.frame(dat[,1:4], umap, tsne$Y, as.factor(g$membership) )
colnames(df) = c(paste0("PC_",1:3), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
df$hotspot <- rownames(dat)


p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("11_whole_genome_leiden_pca",pcn,"_k",k,"_resolution",resolution,"_tsne.png"), p, width = 6, height = 5.3)

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_d3("category20") + 
  labs(color="Cluster")+
  # guides(fill = guide_legend(title = "Cluster"))+
  theme_bw()
ggsave(paste0("11_whole_genome_leiden_pca",pcn,"_k",k,"_resolution",resolution,"_umap.png"), p, width = 3, height = 2.1)

dup <-data.frame(dat_dup[,1:3],hotspot=rownames(dat_dup))
colnames(dup)[1:3] <-paste0("PC_",1:3)
dup <-left_join(dup,df[,1:8],by=c("PC_1","PC_2","PC_3"))
dup <-dup[,c(1:3,5:9,4)]
fdatl <- bind_rows(dup,df)

write.table(fdatl,paste0("11_whole_genome_leiden_pca",pcn,"_k",k,"_resolution",resolution,".txt"),row.names = F, col.names = T,quote =F,sep="\t")
