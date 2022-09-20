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
library(umap)
library(uwot)
library(kernlab)
setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/")
org<-read.csv("6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure")
# load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/emplambda_0/06_permutation_test_1000_sig_kmer.Rdata")
load("../06_permutation_test_1000_sig_kmer.Rdata")
sigK <-fdat0

colnames(org)[1] <-"hotspot"
org$hotspot <-gsub(">","",org$hotspot)
rownames(org)=org$hotspot
org <-org[,-1]

Sorg <-org[,which(colnames(org) %in% sigK$seq)]
save(Sorg,file ="08_permutation_wilocx_overlap_sig_kmer_0_1000.Rdata")
Sorg_all <-Sorg

rn=88
set.seed(rn)
Sorg <-Sorg_all[sample(1:nrow(Sorg_all), 5000, replace = FALSE),]
random5000 <- Sorg
# save(random5000,file="08_random_5000_in_permutation_wilocx_overlap_sig.Rdata")


setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/sampling/")
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 
dat1 <-t(Sorg)
pca=prcomp(dat1)
library(factoextra)
p1 <-fviz_screeplot(pca,ncp=30,addlabels = TRUE,font.x=c(5),font.y=c(16),font.main=c(5))
ggsave(paste0("08_PCA_screeplot_",rn,".pdf"),p1,dpi=300,width=10,height=5)

pca_result = pca$rotation
df <-data.frame(pca_result) 
p =ggplot(df) + geom_point(aes(x = PC1, y = PC2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave(paste0("08_pca_cluster_",rn,".png"),p,dpi=150,width=5,height=5)

tsne = Rtsne::Rtsne(pca_result[,1:30], pca = FALSE,check_duplicates = FALSE)
# save(tsne,file="08_new_pca_tsne.Rdata")

df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave(paste0("08_PCA_tsn_cluster_",rn,".png"),p,dpi=150,width=5,height=5)
#---------------
rumap <- uwot::umap(pca_result[,1:30])
rownames(rumap)=rownames(Sorg)
df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2])
# save(rumap,file="08_new_pca_umap.Rdata")

p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2), size = 0.3, alpha = 0.3)+theme_bw() + p_theme
ggsave(paste0("08_PCA_Umap_cluster_",rn,".png"),p,dpi=150,width=5,height=5)
dat <-Sorg
#-----------------------------------------------------kpca-------------------------------------
library(kernlab)
pca_list <-list()
pca_list[[1]]<- kpca(~.,data=dat, kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)
pca_list[[2]]<- kpca(~.,data=dat, kernel="laplacedot", kpar = list(sigma = 0.01))
pca_list[[3]]<- kpca(~.,data=dat, kernel="tanhdot", kpar = list(scale = 0.01, offset = 1))
pca_list[[4]]<- kpca(~.,data=dat, kernel="besseldot", kpar = list(sigma = 1, order = 1, degree = 1))
pca_list[[5]]<- kpca(~.,data=dat, kernel="anovadot", kpar = list(degree = .1))
pca_list[[6]]<- kpca(~.,data=dat, kernel="vanilladot", kpar = list())
pca_list[[7]]<- kpca(~.,data=dat, kernel="polydot", kpar = list(degree = 2, scale = 1, offset = 5))
pca_list[[8]] <-kpca(~.,data=dat, kernel='splinedot', kpar = list())

names(pca_list)=c("rbfdot","laplacedot","tanhdot","besseldot","anovadot","vanilladot","polydot","splinedot")
# save(pca_list,file = "08_kpca.Rdata")


plot <-function(i=NULL){
    aa=as.data.frame(pca_list[[i]]@rotated)
    colnames(aa)[1:2] <- c("PC1","PC2")
    p =ggplot(aa) + geom_point(aes(x = PC1, y = PC2), size = 0.1, alpha = 0.3)+ggtitle(names(pca_list)[i])+theme_bw() + p_theme
    return(p)
}

plist <-lapply(c(1:8),plot)
pdf(paste0("./08_kpca_sampling_result_",rn,".pdf"),width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()

#-------------------------umap 
uumap <-function(i=NULL){
    aa=as.data.frame(pca_list[[i]]@rotated)
    rumap <- uwot::umap(aa, n_neighbors = 30)
    return(rumap)
    print(i)
}

# umap_list<-lapply(c(1:8),uumap)
# names(umap_list)=c("rbfdot","laplacedot","tanhdot","besseldot","anovadot","vanilladot","polydot","splinedot")
# uplot <-function(i=NULL){
#     df <-data.frame(UMAP_1=umap_list[[i]][,1],UMAP_2=umap_list[[i]][,2])
#     p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.3)+ggtitle(names(umap_list)[i])+theme_bw() + p_theme
#     return(p)
# }

umap_list<-lapply(c(1:7),uumap)
names(umap_list)=c("rbfdot","laplacedot","tanhdot","besseldot","anovadot","vanilladot","polydot")
uplot <-function(i=NULL){
    df <-data.frame(UMAP_1=umap_list[[i]][,1],UMAP_2=umap_list[[i]][,2])
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.3)+ggtitle(names(umap_list)[i])+theme_bw() + p_theme
    return(p)
}

plist <-lapply(c(1:7),uplot)
pdf("./08_kpca_umap_sampling_result.pdf",width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()
#--------------------------------------------tsne
tsne = Rtsne::Rtsne(pca_result[,1:30], pca = FALSE,check_duplicates = FALSE)
ttsne <-function(i=NULL){
    aa=as.data.frame(pca_list[[i]]@rotated)
    tsne <- Rtsne::Rtsne(aa[,1:30], pca = FALSE,check_duplicates = FALSE)
    return(tsne)
    print(i)
}
tsne_list<-lapply(c(1:4,6:8),ttsne)
names(tsne_list)=c("rbfdot","laplacedot","tanhdot","besseldot","vanilladot","polydot","splinedot")

tplot <-function(i=NULL){
    df = data.frame(tsne1 = tsne_list[[i]]$Y[,1], tsne2 = tsne_list[[i]]$Y[,2])
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2), size = 0.3, alpha = 0.3)+ggtitle(names(tsne_list)[i])+theme_bw() + p_theme
    return(p)
}
plist <-lapply(c(1:7),tplot)
pdf("./08_kpca_tsne_sampling_result.pdf",width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()


#----------------------
#---------------------------plot 3d pca
library("gg3D")

plot_3D <-function(i=NULL){
  aa=as.data.frame(pca_list[[i]]@rotated)
  colnames(aa)[1:3] <- c("PC1","PC2","PC3")
  p =ggplot(aa,aes(x = PC1, y = PC2,z=PC3)) + axes_3D()+stat_3D(size=0.001, alpha = 0.3)+
    ggtitle(names(pca_list)[i])+
    theme_void() +
    labs_3D(labs=c("x", "y", "z"),angle=c(0,0,0))#+ p_theme
  return(p)
}

plist <-lapply(c(1:8),plot_3D)
pdf("./08_kpca_sampling_3D_result.pdf",width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()

#_----------------------------------------------------------------------------------------------------------------------------




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

save(re,file="08_sampling_kmeans_3_10.Rdata")
j=1
plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df = data.frame(tsne1 = tsne_list[[j]]$Y[,1], tsne2 = tsne_list[[j]]$Y[,2], cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.8, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}


plist <- lapply(c(1:8),plot)
p =gridExtra::marrangeGrob(plist,ncol=4,nrow=2,top=names(tsne_list[j]))
# p =CombinePlots(plist,ncol=4,nrow=2)
ggsave("08_sampling_kmeans_tsne.png",p,dpi=150,width=12,height=4.8)
#_-------------------------------------------------------------------------

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=umap_list[[j]][,1],UMAP_2=umap_list[[j]][,2],cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+2)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:8),plot)
p =gridExtra::marrangeGrob(plist,ncol=4,nrow=2,top=names(tsne_list[j]))
ggsave("08_sampling_kmeans_UMAP.png",p,dpi=150,width=12,height=4.8)

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

df <-data.frame(UMAP_1=umap_list[[j]][,1],UMAP_2=umap_list[[j]][,2],cluster=g$membership)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("08_random_5000_euc_leiden_Umap_cluster_new.png",p,dpi=300,width=5,height=5)
#----------------------------------------------------------------------------------cluster_louvain
g <-cluster_louvain(g1, weights = NULL)
df <-data.frame(UMAP_1=umap_list[[j]][,1],UMAP_2=umap_list[[j]][,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("08_random_5000_euc_louvain_Umap_cluster_new.png",p,dpi=300,width=5,height=5)

df <-data.frame(tsne1 = tsne_list[[j]]$Y[,1], tsne2 = tsne_list[[j]]$Y[,2],cluster=g$membership)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
ggsave("08_random5000__euc_louvain_tsn_cluster_new.png",p,dpi=300,width=5,height=4.5)

#-----------------------------------------------------------------------------------specc
library("kernlab")

sc <-specc(as.matrix(dat[1:100,]),centers=2)
aa  <-specc(as.matrix(dat[1:100,]),centers=2)

set.seed(1)
dat_1000 <- dat[sample(1:nrow(dat), 1000, replace = FALSE),]
pca=kpca(~.,data=dat_1000 , kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)
pca_result = pca@rotated
tsne = Rtsne::Rtsne(pca_result[,1:30], pca = FALSE,check_duplicates = FALSE)
#---------------
rumap <- uwot::umap(pca_result[,1:30])
rownames(rumap)=rownames(dat_1000)


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