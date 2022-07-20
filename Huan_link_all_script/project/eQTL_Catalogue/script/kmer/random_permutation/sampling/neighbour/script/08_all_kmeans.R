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

setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/")
# org<-read.csv("6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
# setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/")
# load("06_wilcox_permutation_overlap_significant.Rdata")
# sigK <-wilcox_permutation_overlap

# colnames(org)[1] <-"hotspot"
# org$hotspot <-gsub(">"," ",org$hotspot)
# rownames(org)=org$hotspot
# org <-org[,-1]


setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/")
# Sorg <-org[,which(colnames(org) %in% sigK$seq)]
# save(Sorg,file ="08_permutation_wilocx_overlap_sig_kmer_0_100.Rdata")
load("08_permutation_wilocx_overlap_sig_kmer_0_100.Rdata")

Sorg_all <-Sorg



setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/")
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 
dat1 <-t(Sorg)
pca=prcomp(dat1)
save(pca,file="08_all_PCA.Rdata")
library(factoextra)
p1 <-fviz_screeplot(pca,ncp=30,addlabels = TRUE,font.x=c(5),font.y=c(16),font.main=c(5))
ggsave("08_all_PCA_screeplot.pdf",p1,dpi=300,width=10,height=5)

pca_result = pca$rotation
df <-data.frame(pca_result) 
p =ggplot(df) + geom_point(aes(x = PC1, y = PC2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("08_all_pca_cluster_new.png",p,dpi=150,width=5,height=5)

tsne = Rtsne::Rtsne(pca_result[,1:30], pca = FALSE,check_duplicates = FALSE)
save(tsne,file="08_all_new_pca_tsne.Rdata")

df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("08_all_PCA_tsn_cluster_new.png",p,dpi=150,width=5,height=5)
#---------------
rumap <- uwot::umap(pca_result[,1:30])
rownames(rumap)=rownames(Sorg)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2])
save(rumap,file="08_all_new_pca_umap.Rdata")

p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("08_all_PCA_Umap_cluster_new.png",p,dpi=150,width=5,height=5)
dat <-Sorg

#------------------------------------------------------------------------------kmeans
library(amap)
dat1 <-as.matrix(dat)
kr_n <-function(i){
    print(paste(i,"start",sep="\t"))
    kr <-Kmeans(dat1, i)
    print(paste(i,"end",sep="\t"))
    return(kr)
}

re <-lapply(c(4:6),kr_n)

save(re,file="08_all_kmeans_4_6.Rdata")
plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.8, alpha = 0.3)+ggtitle(i+3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}


plist <- lapply(c(1:3),plot)
p =CombinePlots(plist,ncol=3,nrow=1)
ggsave("08_all_kmeans_tsne.png",p,dpi=150,width=12,height=3.5)
#_-------------------------------------------------------------------------

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:3),plot)
p =CombinePlots(plist,ncol=3,nrow=1)
ggsave("08_all_kmeans_UMAP.png",p,dpi=150,width=12,height=3.5)

#-------------------------------------------------------------------------------------------------leiden
library(parallelDist)
# eu = parallelDist::parDist(x = as.matrix(dat), method = "euclidean",threads = 20)
# set.seed(1234)
# sample_d <- sample(1:nrow(dat), 40000, replace = FALSE)
# sample_dat <-dat[sample_d,]
sample_dat <-dat
eu = parallelDist::parDist(x = as.matrix(sample_dat), method = "euclidean",threads = 20)
# set.seed(1)
# r_a <-sample(1:length(eu), 20000000, replace = FALSE)
# eu_d <-eu[r_a]
# eu_df <-data.frame(eu_d)
cutoff <-quantile(eu,0.75) #22.17723
names(cutoff)=NULL

tu <-eu
tu[which(tu<=cutoff)]=0
tu[which(tu>cutoff)]=1
library(igraph)
g1 <- graph_from_adjacency_matrix(tu,mode="undirected")

g <-cluster_leiden(g1, resolution_parameter=0.03)
leiden_r <-data.frame(hotspot=rownames(sample_dat),cluster=g$membership)

df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("08_random_5000_euc_leiden_Umap_cluster_new.png",p,dpi=300,width=5,height=5)
#----------------------------------------------------------------------------------cluster_louvain
g <-cluster_louvain(g1, weights = NULL)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("08_random_5000_euc_louvain_Umap_cluster_new.png",p,dpi=300,width=5,height=5)

df <-data.frame(tsne1=tsne$Y[,1],tsne2=tsne$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("08_random5000__euc_louvain_tsn_cluster_new.png",p,dpi=300,width=5,height=4.5)

#-----------------------------------------------------------------------------------specc
library("kernlab")

sc <-specc(as.matrix(dat[1:100,]),centers=2)
aa  <-specc(as.matrix(dat[1:100,]),centers=2)

set.seed(1)
dat_1000 <- dat[sample(1:nrow(dat), 1000, replace = FALSE),]



specc_n <-function(i){
    print(paste(i,"start",sep="\t"))
    sc <-specc(as.matrix(dat_1000),centers=i)
    print(paste(i,"end",sep="\t"))
    return(sc)
}
re <-lapply(c(3:10),specc_n)
# re <-mclapply(c(3:10),specc_n,mc.cores=9)
#---------------
# dat1 <-t(Sorg)
pca=prcomp(t(dat_1000))
pca_result = pca$rotation
tsne = Rtsne::Rtsne(pca_result[,1:30], pca = FALSE,check_duplicates = FALSE)
#---------------
rumap <- uwot::umap(pca_result[,1:30])
rownames(rumap)=rownames(dat_1000)


plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]@.Data)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("08_random_1000_specc_UMAP.png",p,dpi=150,width=12,height=5)

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2],cluster = re[[i]]@.Data)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}
plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("08_random_1000_specc_tsne.png",p,dpi=150,width=12,height=5)


#------------------------------------------------------------------------------------------------------------mcluster

mc <-function(i){
    print(paste(i,"start",sep="\t"))
    mod <-Mclust(dat_1000,G=1:i)
    print(paste(i,"end",sep="\t"))
    return(mod)
}
models <-lapply(c(3:10),mc)
plot <-function(i){
    re <-data.frame(models[[i]]$classification)
    colnames(re)[1]="cluster"
    df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.8, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    return(p)
}
plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("08_random1000_mclust_tsne_new.png",p,dpi=120,width=12,height=5)

uplot <-function(i=NULL){
    re <-data.frame(models[[i]]$classification)
    colnames(re)[1]="cluster"
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
}

plist <- lapply(c(1:8),uplot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("08_random1000_mclust_UMAP_new.png",p,dpi=120,width=12,height=5)