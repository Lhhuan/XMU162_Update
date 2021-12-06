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
library(parallel)
library(pheatmap)
setwd("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_6/kmer/6/")

org<-read.csv("6mers.csv",header = T,sep = ",") %>% as.data.frame()
colnames(org)[1] <-"hotspot"
org2 <-melt(org,"hotspot")
colnames(org2)[2] <-"seq"
org2$hotspot <-str_replace(org2$hotspot,">","")

cluster_result <-read.table("communities.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(cluster_result)<-c("chr","start","end","community")
cluster_result$hotspot <-paste0(cluster_result$chr,":",cluster_result$start,"-",cluster_result$end)

# rownames(cluster_result)<-cluster_result$hotspot
cluster_result<-cluster_result%>%select(hotspot,community)

m <- dcast(data=org2,hotspot~seq)
aaa <-inner_join(m,cluster_result,by="hotspot")
bbb <-aaa[order(aaa$community),]
final_fig_used <- bbb%>%select(-community)
rownames(final_fig_used) <-final_fig_used[,1]
n <-final_fig_used[,-1]
annotation_rows <- bbb%>%select(community,hotspot)
rownames(annotation_rows)<-annotation_rows$hotspot
annotation_rows<-annotation_rows%>%select(-hotspot)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
f_n<-apply(n,2,range01)%>%data.frame() 

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_6/figure/")
pheatmap(f_n,cluster_col = FALSE, cluster_rows = FALSE,border=FALSE,
         annotation_row = annotation_rows,show_colnames = F,show_rownames = F,fontsize_row=5.5,
         annotation_legend = T,annotation_names_row = T, annotation_names_col = F,scale="none",
         color = colorRampPalette(c("white","#3edbf0"))(50),
         file='chr1_6_6kmer_5_community.png')

print("finish")

# pheatmap(n[1:300,1:500],cluster_col = FALSE, cluster_rows = FALSE,border=FALSE,
#          annotation_row = annotation_rows,show_colnames = F,show_rownames = F,
#          annotation_legend = F,annotation_names_row = F, annotation_names_col = F,scale="column",
#          file='chr1_6_6kmer_5_community.png')

# pheatmap(n[1:300,1:300],file="test_row.png",cluster_col = FALSE, cluster_rows = FALSE,border=FALSE,show_colnames = F,show_rownames = F,scale="row")


# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # n_0_1<-range01(n[1:300,1:300])
# n_0_1<-apply(n[1:300,1:300],2,range01)%>%data.frame()  
# pheatmap(n_0_1,file="test_none.png",cluster_col = FALSE, cluster_rows = FALSE,border=FALSE,show_colnames = F,show_rownames = F,scale="none")


# rownames(cluster_result)<-cluster_result$hotspot
# cluster_result<-cluster_result%>%select(community)