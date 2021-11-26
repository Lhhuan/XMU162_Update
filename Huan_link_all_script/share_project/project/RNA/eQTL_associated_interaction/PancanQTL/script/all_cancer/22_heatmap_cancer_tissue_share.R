library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(pheatmap)
library(reshape2)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/cancer_total/")

org <- read.table("21_cancer_pair_overlap_index.txt.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(org) <-c("Index","tissue1","tissue2")
org2 <-org[,c(1,3,2)]
colnames(org2)[2] ="tissue1"
colnames(org2)[3] ="tissue2"
org_final <-unique(bind_rows(org,org2))
org_final =org_final[,c(2,3,1)]
m <- dcast(data=org_final,tissue1~tissue2)
rownames(m) <-m[,1]
n <-m[,-1]
save(n,file="cancer_cancer_int.Rdata")
names_tissue <-colnames(n)
 
ann_colors = list(
  tissue = c(ACC = "#61b15a",BRCA="#28abb9", COAD="#ffab73", ESCA="#c19277", KICH="#ccffbd", KIRC="#ccffbd", KIRP="#ccffbd", LIHC="#beca5c", LUAD="#c0e218",LUSC="#c0e218", OV="#ffc1fa", PAAD="#87431d",PRAD="#e8eae6", SKCM="#845ec2", STAD="#eac8af", TGCT="#9088d4", THCA="#064420", UCEC="#f09ae9",UCS="#f09ae9", LAML="#de8971")   
)



annotation_row = data.frame(tissue = factor(colnames(n)))
rownames(annotation_row) <- annotation_row$tissue
annotation_col = data.frame(tissue = factor(colnames(n)))
rownames(annotation_col) <- annotation_col$tissue
annotation_row$tissue <- gsub("-", "_", annotation_row$tissue)
annotation_col$tissue <- gsub("-", "_", annotation_col$tissue)
# pheatmap(n,cluster_col = FALSE, cluster_rows = FALSE,border=FALSE, annotation_col = annotation_col,
#          annotation_row = annotation_row,annotation_colors = ann_colors,
#          cellwidth = 10, cellheight = 6, main = "NAR")
# #------------
# pdf("sssss.pdf",width = 10,height=10)
pheatmap(n,cluster_col = FALSE, cluster_rows = FALSE,border=FALSE, annotation_col = annotation_col,
         annotation_row = annotation_row,annotation_colors = ann_colors,show_colnames = F,fontsize_row=5.5,
         annotation_legend = F,annotation_names_row = F, annotation_names_col = F,
         cellwidth = 12, cellheight = 12, scale="none",
         color = colorRampPalette(c("#3edbf0","white",  "#ff4646"))(50),
         file='22_heatmap_cancer_cancer_share_hotspot_before.pdf')

pheatmap(n,cluster_col = T, cluster_rows = T,border=FALSE, annotation_col = annotation_col,
         annotation_row = annotation_row,annotation_colors = ann_colors,show_colnames = F,fontsize_row=5.5,
         annotation_legend = F,annotation_names_row = F, annotation_names_col = F,
         cellwidth = 12, cellheight = 12, scale="none",
         color = colorRampPalette(c("#cceeff","#33bbff"))(50),
         file='22_heatmap_cancer_cancer_share_hotspot_cluster.pdf')

pheatmap(n,cluster_col = FALSE, cluster_rows = FALSE,border=FALSE, annotation_col = annotation_col,
         annotation_row = annotation_row,annotation_colors = ann_colors,show_colnames = F,fontsize_row=5.5,
         annotation_legend = F,annotation_names_row = F, annotation_names_col = F,
         cellwidth = 12, cellheight = 12, scale="none",
         color = colorRampPalette(c("#cceeff","#33bbff"))(50),
         file='22_heatmap_cancer_cancer_share_hotspot.pdf')


# library(pheatmap)
# # a<-read.table("COO.txt",sep="\t",header=T)
# # a<-a[order(a$Predicted.classification),]
# # table(a$Predicted.classification)
# # annotation_col = data.frame(Cell_Type=(c('ABC','GCB')[a$Predicted.classification]))
# # rownames(annotation_col) = paste("", 1:441, sep = "")
# # pheatmap(b,treeheight_row=100,treeheight_col=20,cluster_cols=TRUE,color=colorRampPalette(c("green","black","red"))(1000),
# #  border_color=NA,fontsize=10,fontsize_row=8,fontsize_col=16,file='111.jpg',annotation_col = annotation_col, scale = "row")
# annotation_col = data.frame(CellType = factor(c(rep("exp_D3",3), rep("exp_H", 3)))
# rownames(annotation_col) = paste("Test", 1:10, sep = "")
# head(annotation_col)
# heatmap <- pheatmap(a,scale = "row",cluster_col = FALSE, border=FALSE, treeheight_row = 30,cellwidth = 30, cellheight = 12, main = "Glycolysis PMID:24140020")