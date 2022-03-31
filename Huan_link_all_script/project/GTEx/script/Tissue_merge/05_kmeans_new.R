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
library(Rtsne)
library(irlba)

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

# save(dat,file="04_hotspot_mark_fraction.Rdata")
library(amap)
dat1 <-as.matrix(dat)
kr_n <-function(i){
    print(paste(i,"start",sep="\t"))
    kr <-Kmeans(dat1, i)
    print(paste(i,"end",sep="\t"))
    return(kr)
}

re <-lapply(c(2:10),kr_n)


load("04_new_pca_tsne.Rdata")
load("04_new_pca_umap.Rdata")
# tsne = Rtsne::Rtsne(pca_result, pca = FALSE,check_duplicates = FALSE)

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.8, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}


plist <- lapply(c(1:9),plot)
p =CombinePlots(plist,ncol=3,nrow=3)
ggsave("05_kmeans_tsne_new.png",p,dpi=300,width=11.5,height=10)
#_-------------------------------------------------------------------------

# rumap <- uwot::umap(pca_result, n_neighbors = 30)

plot <-function(i){
    print(paste(i,"start",sep="\t"))
    df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = re[[i]]$cluster)
    p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle(i+1)+theme_bw() + p_theme+
    guides(col=guide_legend(title='Cluster'))
    print(paste(i,"end",sep="\t"))
    return(p)
}

plist <- lapply(c(1:9),plot)
p =CombinePlots(plist,ncol=3,nrow=3)
ggsave("05_kmeans_UMAP_new.png",p,dpi=300,width=11.5,height=10)


library(clusterCrit)
intCriteria(dat1, kr$cluster, 'Calinski_Harabasz')
#24514.25