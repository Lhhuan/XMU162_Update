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

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/")
org<-read.csv("6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()

sigK <-read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/sampling/emplambda_0/figure/hotspot_emplambda_0_random_kmer_wilcox.txt",header = T,sep = "\t") %>% as.data.frame()
sigK <-filter(sigK,annotation !="ns")
colnames(org)[1] <-"hotspot"
org$hotspot <-gsub(">"," ",org$hotspot)
rownames(org)=org$hotspot
org <-org[,-1]

Sorg <-org[,which(colnames(org) %in% sigK$seq)]
setwd("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/figure/emp0_100/")
save(Sorg,file ="sig_kmer_0_100.Rdata")

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 
dat1 <-t(Sorg)
pca=prcomp(dat1)
library(factoextra)
p1 <-fviz_screeplot(pca,ncp=30,addlabels = TRUE,font.x=c(5),font.y=c(16),font.main=c(5))
ggsave("08_PCA_screeplot.pdf",p1,dpi=300,width=10,height=5)

pca_result = pca$rotation
df <-data.frame(pca_result) 
p =ggplot(df) + geom_point(aes(x = PC1, y = PC2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("08_pca_cluster_new.png",p,dpi=150,width=5,height=5)

tsne = Rtsne::Rtsne(pca_result[,1:30], pca = FALSE,check_duplicates = FALSE)
save(tsne,file="08_new_pca_tsne.Rdata")

df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("08_PCA_tsn_cluster_new.png",p,dpi=150,width=5,height=5)
#---------------
rumap <- uwot::umap(pca_result[,1:30])
rownames(rumap)=rownames(Sorg)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2])
save(rumap,file="08_new_pca_umap.Rdata")

p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave("08_PCA_Umap_cluster_new.png",p,dpi=150,width=5,height=5)
dat <-Sorg
#-----------mcluster


mc <-function(i){
    print(paste(i,"start",sep="\t"))
    mod <-Mclust(dat,G=1:i)
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
ggsave("08_mclust_tsne_new.png",p,dpi=120,width=12,height=6)

uplot <-function(i=NULL){
    re <-data.frame(models[[i]]$classification)
    colnames(re)[1]="cluster"
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
}

plist <- lapply(c(1:8),uplot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("08_mclust_UMAP_new.png",p,dpi=120,width=12,height=6)


#------------------------------------------------------------------------------kmeans
library(amap)
dat1 <-as.matrix(dat)
kr_n <-function(i){
    print(paste(i,"start",sep="\t"))
    kr <-Kmeans(dat1, i)
    print(paste(i,"end",sep="\t"))
    return(kr)
}

re <-lapply(c(3:10),kr_n)

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.8, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}


plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("08_kmeans_tsne.png",p,dpi=150,width=12,height=6)
#_-------------------------------------------------------------------------

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("05_kmeans_UMAP.png",p,dpi=150,width=12,height=6)

#-------------------------------------------------------------------------------------------------leiden
library(parallelDist)
eu = parallelDist::parDist(x = as.matrix(dat), method = "euclidean",threads = 20)
set.seed(1)
r_a <-sample(1:length(eu), 20000000, replace = FALSE)
eu_d <-eu[r_a]
eu_df <-data.frame(eu_d)
cutoff <-quantile(eu_d,0.75) #76.51333
names(cutoff)=NULL

set.seed(111)
test <- sample(1:nrow(dat),size=10000,replace =T)
dat_test <-dat[test,]
tu <-parallelDist::parDist(x = as.matrix(dat_test), method = "euclidean",threads = 10)
tu[which(tu<=cutoff)]=0
tu[which(tu>cutoff)]=1
library(igraph)
g1 <- graph_from_adjacency_matrix(tu,mode="undirected")

g <-cluster_leiden(g1, resolution_parameter=0.03)
leiden_r <-data.frame(hotspot=rownames(dat_test),cluster=g$membership)
dat1 <-t(dat_test)
pca=prcomp(dat1)
pca_result = pca$rotation

rumap <- uwot::umap(pca_result[,1:30])
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("08_leiden_Umap_cluster_new.png",p,dpi=300,width=5,height=5)
#----------------------------------------------------------------------------------cluster_louvain
g <-cluster_louvain(g1, weights = NULL)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("08_louvain_Umap_cluster_new.png",p,dpi=300,width=5,height=5)


#-----------------------------------------------------------------------------------specc
library("kernlab")

sc <-specc(as.matrix(dat[1:100,]),centers=2)
aa  <-specc(as.matrix(dat[1:100,]),centers=2)

specc_n <-function(i){
    print(paste(i,"start",sep="\t"))
    sc <-specc(as.matrix(dat),centers=i)
    print(paste(i,"end",sep="\t"))
}
re <-lapply(c(3:10),specc_n)
# re <-mclapply(c(3:10),specc_n,mc.cores=9)
plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:8),plot)
p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("05_kmeans_UMAP.png",p,dpi=150,width=12,height=6)

write.table(hotspot_kmer_need_test_value,"./figure/hit_hospot_ratio_kmer_value.txt",row.names = F, col.names = T,quote =F,sep="\t")

gzip("./figure/hit_hospot_ratio_kmer_value.txt")