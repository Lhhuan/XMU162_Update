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
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) #

#---------------------------------------------------------
library(Hmisc)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")

org<-read.table("../../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$chr,org$start,org$end,sep="_")
org1 <-org%>%group_by(Marker,hotspot)%>%summarise(sum_overlap_fraction =sum(overlap_fraction))%>%data.frame()
dat  <-dcast(org1[,c(2,1,3)],hotspot~Marker)

rownames(dat) = dat[,"hotspot"]
dat = dat[,-1]

# mod1 <- Mclust(dat)
# re <-data.frame(mod1$classification)
# tsne = Rtsne::Rtsne(dat, pca = FALSE,check_duplicates = FALSE)

# df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re$mod1.classification)
# p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.8, alpha = 0.3)+ggtitle("Mclust")+theme_bw() + p_theme+
# guides(col=guide_legend(title='Cluster'))
# ggsave("05_mclust.png",p,dpi=150,width=6,height=5)

# mod2 <- Mclust(dat,G=1:7)
mc <-function(i){
    print(paste(i,"start",sep="\t"))
    mod <-Mclust(dat,G=1:i)
    print(paste(i,"end",sep="\t"))
    return(mod)
}

models <-lapply(c(2:10),mc)

load("04_new_pca_tsne.Rdata")
load("04_new_pca_umap.Rdata")

plot <-function(i){
    re <-data.frame(models[[i]]$classification)
    colnames(re)[1]="cluster"
    df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.8, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    return(p)
}
plist <- lapply(c(1:9),plot)
p =CombinePlots(plist,ncol=3,nrow=3)
ggsave("05_mclust_tsne_new.png",p,dpi=120,width=10,height=10)

#-------------------------------------------------------------

uplot <-function(i=NULL){
    re <-data.frame(models[[i]]$classification)
    colnames(re)[1]="cluster"
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
}

plist <- lapply(c(1:9),uplot)
p =CombinePlots(plist,ncol=3,nrow=3)
ggsave("05_mclust_UMAP_new.png",p,dpi=120,width=10,height=10)


# df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = result$cluster_df$cluster)
# p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle("KKL")+theme_bw() + p_theme
# ggsave("05_KKL_Umap_cluster.png",p,dpi=150,width=5.5,height=5)

