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
library(umap)
library(uwot)
library(kernlab)
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 

#---------------------------------------------------------
library(Hmisc)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")

org<-read.table("../../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$chr,org$start,org$end,sep="_")
org <-org[,c("Marker","overlap_fraction","hotspot")]
org1 <-org%>%group_by(Marker,hotspot)%>%summarise(sum_overlap_fraction =sum(overlap_fraction))%>%data.frame()
# org1%>%group_by(Marker)%>%summarise(count=n())
org2 <-dcast(org1[,c(2,1,3)],hotspot~Marker)
rownames(org2)<-org2$hotspot
dat <-org2[,-1]
set.seed(111)
test <- sample(1:nrow(dat),size=10000,replace =T)
pca_list <-list()
pca_list[[1]]<- kpca(~.,data=dat[test,], kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)
pca_list[[2]]<- kpca(~.,data=dat[test,], kernel="laplacedot", kpar = list(sigma = 0.01))
pca_list[[3]]<- kpca(~.,data=dat[test,], kernel="tanhdot", kpar = list(scale = 0.01, offset = 1))
pca_list[[4]]<- kpca(~.,data=dat[test,], kernel="besseldot", kpar = list(sigma = 1, order = 1, degree = 1))
pca_list[[5]]<- kpca(~.,data=dat[test,], kernel="anovadot", kpar = list(degree = .1))
pca_list[[6]]<- kpca(~.,data=dat[test,], kernel="vanilladot", kpar = list())
pca_list[[7]]<- kpca(~.,data=dat[test,], kernel="polydot", kpar = list(degree = 2, scale = 1, offset = 5))
pca_list[[8]] <-kpca(~.,data=dat[test,], kernel='splinedot', kpar = list())

names(pca_list)=c("rbfdot","laplacedot","tanhdot","besseldot","anovadot","vanilladot","polydot","splinedot")
save(pca_list,file = "04_kpca.Rdata")


plot <-function(i=NULL){
    aa=as.data.frame(pca_list[[i]]@rotated)
    colnames(aa)[1:2] <- c("PC1","PC2")
    p =ggplot(aa) + geom_point(aes(x = PC1, y = PC2), size = 0.1, alpha = 0.3)+ggtitle(names(pca_list)[i])+theme_bw() + p_theme
    return(p)

}

plist <-lapply(c(1:8),plot)
pdf("./04_kpca_sampling_result.pdf",width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()

#-------------------------umap 
uumap <-function(i=NULL){
    aa=as.data.frame(pca_list[[i]]@rotated)
    rumap <- uwot::umap(aa, n_neighbors = 30)
    return(rumap)
    print(i)
}

umap_list<-lapply(c(1:8),uumap)
names(umap_list)=c("rbfdot","laplacedot","tanhdot","besseldot","anovadot","vanilladot","polydot","splinedot")
uplot <-function(i=NULL){
    df <-data.frame(UMAP_1=umap_list[[i]][,1],UMAP_2=umap_list[[i]][,2])
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.3)+ggtitle(names(umap_list)[i])+theme_bw() + p_theme
    return(p)
}

plist <-lapply(c(1:8),uplot)
pdf("./04_kpca_umap_sampling_result.pdf",width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()

#---------------------------plot 3d pca
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
pdf("./04_kpca_sampling_3D_result.pdf",width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()



#------------------











# #----------------------rbfdot
# kpc <-kpca(~.,data=dat[test,],kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)
# pdf("05_rbfdot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()
# #----------------laplacedot
# kpc <-kpca(~.,data=dat[test,],kernel="laplacedot", kpar = list(sigma = 0.01))
# pdf("05_laplacedot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()
# #----------------tanhdot
# kpc <- kpca(~.,data=dat[test,], kernel="tanhdot", kpar = list(scale = 0.01, offset = 1))

# pdf("05_tanhdot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()
# #-------------------------besseldot
# kpc <- kpca(~.,data=dat[test,], kernel="besseldot", kpar = list(sigma = 1, order = 1, degree = 1))
# pdf("05_besseldot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()
# #--------------------------anovadot
# kpc <- kpca(~.,data=dat[test,], kernel="anovadot", kpar = list(degree = .1))
# pdf("05_anovadot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()
# #---------------------vanilladot liner
# kpc <- kpca(~.,data=dat[test,], kernel="vanilladot", kpar = list())
# pdf("05_vanilladot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()
# #-----------------------polydot
# kpc <-kpca(~.,data=dat[test,], kernel="polydot", kpar = list(degree = 2, scale = 1, offset = 5))
# pdf("05_polydot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()
# #----------------------
# #----------------------splinedot
# kpc <-kpca(~.,data=dat[test,], kernel='splinedot', kpar = list())
# pdf("05_splinedot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()
# #----------------
# kpc <-kpca(~.,data=dat[test,],kernel = "stringdot", kpar = list( lambda = 0.5))

# pdf("05_stringdot.pdf")
# plot(rotated(kpc),
#      xlab="1st Principal Component",ylab="2nd Principal Component")
# dev.off()


# Gpca <-kpca(as.matrix(dat[test,]),kernel="rbfdot",kpar = list(sigma = 0.1),eatures = 0, th = 1e-4)



# GC<-read.table("../../../output/Tissue_merge/Cis_eQTL/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge_filter_3103833_GC.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# colnames(GC)[1:5] <-c("chr","start","end","AT_content","GC_content")
# GC$hotspot <-paste(GC$chr,GC$start,GC$end,sep="_")
# GC <-GC[,c("hotspot","GC_content")]

# dat <- inner_join(GC,org2,by="hotspot")
# rownames(dat) = dat[,"hotspot"]
# dat = dat[,-1]

# # save(dat,file="04_hotspot_mark_fraction.Rdata")
# library(irlba)
# set.seed(11)
# pca = prcomp_irlba(dat, n = 10)
# pca_result = pca$x
# rownames(pca_result) = rownames(dat)
# df <-data.frame(pca_result) 
# p =ggplot(df) + geom_point(aes(x = PC1, y = PC2), size = 0.3, alpha = 0.3)+ggtitle("Add GC PCA")+theme_bw() + p_theme
# ggsave("04_pca_cluster_GC.png",p,dpi=150,width=5,height=5)

# tsne = Rtsne::Rtsne(pca$x, pca = FALSE,check_duplicates = FALSE)

# df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])
# p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2), size = 0.3, alpha = 0.3)+ggtitle("Add GC PCA")+theme_bw() + p_theme
# ggsave("04_tsn_cluster_GC.png",p,dpi=150,width=5,height=5)




# # df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re$mod1.classification)
# # p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.4, alpha = 0.3)+ggtitle("Mclust")+theme_bw() + p_theme+
# # guides(col=guide_legend(title='Cluster'))
# # ggsave("05_mclust.png",p,dpi=150,width=6,height=5)

# #---------------


# dat=dat[,-1] #-GC_content

# pca = prcomp_irlba(dat, n = 9)
# pca_result = pca$x
# rownames(pca_result) = rownames(dat)
# df <-data.frame(pca_result) 
# p =ggplot(df) + geom_point(aes(x = PC1, y = PC2), size = 0.3, alpha = 0.3)+ggtitle("Add GC PCA")+theme_bw() + p_theme
# ggsave("04_pca_cluster.png",p,dpi=150,width=5,height=5)

# tsne = Rtsne::Rtsne(pca$x, pca = FALSE,check_duplicates = FALSE)

# df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2])
# p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2), size = 0.3, alpha = 0.3)+ggtitle("10 factors")+theme_bw() + p_theme
# ggsave("04_tsn_cluster.png",p,dpi=150,width=5,height=5)
# #---------------
# rumap <- uwot::umap(pca_result, n_neighbors = 30)
# df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2])

# p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2), size = 0.3, alpha = 0.3)+ggtitle("10 factors")+theme_bw() + p_theme
# ggsave("04_Umap_cluster.png",p,dpi=150,width=5,height=5)
# #------------------------------------

# eu_dist <- dist(dat, method = "euclidean",diag = TRUE,upper = T)
# # g <-graph_from_adjacency_matrix(eu_dist)
# # df <- melt(as.matrix(eu_dist), varnames = c("hotspot1", "hotspot2"))
# # colnames(df)[3] <-"euclidean_dist"
# # save(df,file="04_10_markers_euclidean_dist_df.Rdata")
# save(eu_dist,file="04_10_markers_euclidean_dist_matrix.Rdata")

