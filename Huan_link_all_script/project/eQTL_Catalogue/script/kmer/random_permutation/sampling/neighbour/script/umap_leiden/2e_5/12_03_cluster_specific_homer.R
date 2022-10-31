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
library(ggsci)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/re_2e_5/homer_200/")

org <-read.table("../../../11_whole_genome_umap_leiden_pca5_k50_resolution2e-05.txt",header = T,sep = "\t") %>% as.data.frame()

merge_motif <-function(i=NULL){
    dat<-read.table(paste0(i,"/knownResults.txt"),header = F,sep = "\t",skip=1) %>% as.data.frame()
    colnames(dat)[1:5] <- c("Motif","Consensus","P_value","Log_P_value","q_value")
    dat <- filter(dat,q_value<0.05)
    if(nrow(dat)>0){
        dat$Motif <-gsub("/.*","",dat$Motif)
        dat$Cluster <-i 
        dat <-dat[,c("Motif","Cluster")]
        return(dat)
    }
    # print(i)
}

tmp <-lapply(unique(org$cluster),merge_motif)
result <-do.call(rbind,tmp)

cluster_n <-result%>%group_by(Cluster)%>%summarize(Number_of_motif=n())%>%data.frame()
cluster_n$Cluster <-as.factor(cluster_n$Cluster)

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

p1 <-ggplot(data = cluster_n, mapping = aes(x =Cluster, y = Number_of_motif,fill=Cluster)) +
     geom_bar(stat = 'identity', width=0.8)+
     scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
     ggtitle("Number of motif")+
     theme(legend.position ="none")+
     labs(y="Number of motif")+
     p_theme
pdf("12_barpolt_number_of_motif.pdf",height=4.7,width=4.5)
print(p1)
dev.off()

cluster_specific <- function(i=NULL){
    cluster <-filter(result,Cluster==i)
    other <-filter(result,Cluster!=i)
    cluster_spe <-filter(cluster,!(Motif%in%other$Motif))
    if(nrow(cluster_spe)>0){
        return(cluster_spe)
    }
}
tmp2 <-lapply(unique(org$cluster),cluster_specific)
result2 <-do.call(rbind,tmp2)
cluster_specifc_n <-result2%>%group_by(Cluster)%>%summarize(Number_of_motif=n())%>%data.frame()

cluster_specifc_n$Cluster <-as.factor(cluster_specifc_n$Cluster)
p1 <-ggplot(data = cluster_specifc_n, mapping = aes(x =Cluster, y = Number_of_motif,fill=Cluster)) +
     geom_bar(stat = 'identity', width=0.8)+
     scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#7F7F7F","#BCBC22","#15BECE"))+
     ggtitle("Number of cluster specific motif")+
     theme(legend.position ="none")+
     labs(y="Number of motif")+
     p_theme
pdf("12_barpolt_number_of_cluster_specific_motif.pdf",height=4.7,width=4.5)
print(p1)
dev.off()

result2<-result2[order(result2$Cluster),]
write.table(result2,"12_cluster_specific_motif.txt",col.names=T,row.names=F,quote=F,sep="\t")

