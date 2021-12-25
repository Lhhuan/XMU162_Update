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
library(mclust)
p_theme<-theme(
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5))

#---------------------------------------------------------
library(Hmisc)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")

org<-read.table("../../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$chr,org$start,org$end,sep="_")
org1 <-org%>%group_by(Marker,hotspot)%>%summarise(sum_overlap_fraction =sum(overlap_fraction))%>%data.frame()
dat  <-dcast(org1[,c(2,1,3)],hotspot~Marker)

rownames(dat) = dat[,"hotspot"]
dat = dat[,-1]

save(dat,file="04_hotspot_10_factor.Rdata")
# eu_dist <- dist(dat, method = "euclidean",diag = TRUE)
# df <-melt(as.matrix(eu_dist))


# a <-as.matrix(eu_dist)
# load("04_10_markers_euclidean_dist_df.Rdata")
# quantile(eu_dist,0.75)

library(parallelDist)
eu = parallelDist::parDist(x = as.matrix(dat), method = "euclidean",threads = 10)
set.seed(1)
r_a <-sample(1:length(eu), 2000000, replace = FALSE)
eu_d <-eu[r_a]
eu_df <-data.frame(eu_d)
# p <-ggplot(eu_df,aes(x=eu_d))+geom_density(color="#5659D4") +theme_bw() #+p_theme
# p
# dev.off()

p <-ggplot(eu_df,aes(x =1, y=eu_d))+geom_violin(color="#5659D4",width=0.3)+geom_boxplot(width=0.2,color="#5659D4")+theme_bw()+p_theme+
labs(x="",y="Dist")+theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())
ggsave("05_dist_cutoff.png",p,width=4,height=4)

quantile(eu_d,0.75)
#1.5933

# adj <-matrix(nrow=nrow(dat),ncol=nrow(dat))
# colnames(adj) <-c(1:122347)
# rownames(adj) <-c(1:122347)
cutoff=1.5933
# eu[which(eu>cutoff)]=1
eu = parallelDist::parDist(x = as.matrix(dat), method = "euclidean",threads = 10)
eu[which(eu<=cutoff)]=0
eu[which(eu>cutoff)]=1

library(igraph)

g1 <- graph_from_adjacency_matrix(eu,mode="undirected")

#-------test leiden
set.seed(111)
test <- sample(1:nrow(dat),size=10000,replace =T)
dat_test <-dat[test,]
tu <-parallelDist::parDist(x = as.matrix(dat_test), method = "euclidean",threads = 10)
tu[which(tu<=cutoff)]=0
tu[which(tu>cutoff)]=1
g1 <- graph_from_adjacency_matrix(tu,mode="undirected")

g <-cluster_leiden(g1, resolution_parameter=0.5)
leiden_r <-data.frame(hotspot=rownames(dat_test),cluster=g$membership)
library(irlba)
pca = prcomp_irlba(dat_test, n = 9)
pca_result = pca$x
rownames(pca_result) = rownames(dat_test)

rumap <- uwot::umap(pca_result, n_neighbors = 30)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("05_leiden_Umap_cluster.png",p,dpi=300,width=5,height=5)
#-----------------------cluster_louvain
g <-cluster_louvain(g1, weights = NULL)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("05_louvain_Umap_cluster.png",p,dpi=300,width=5,height=5)



#------------test
ppp <-head(dat,5)
rownames(ppp) <- c(1:nrow(ppp))


au = parallelDist::parDist(x = as.matrix(ppp), method = "euclidean",threads = 10)
cutoff=0.7
au[which(au>cutoff)]=1
au[which(au<=cutoff)]=0
adt <-matrix(nrow=nrow(ppp),ncol=nrow(ppp))
colnames(adt) <-c(1:5)
rownames(adt) <-c(1:5)
for(i in c(1:nrow(adt))){
    adt[i,i]=0

}


#--------
library(igraph)
adjm <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.9,0.1)), nc=10)
g1 <- graph_from_adjacency_matrix( au)