library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/Knowledge_Graph/output/")

hotspot_g <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(hotspot_g) <-c("chr","start","end","egene")
hotspot_g$Hotspot <-paste(hotspot_g$chr,hotspot_g$start,hotspot_g$end,sep="_")
#------------------nodes
nodes_anno <- read.table("./nodes_annotation/nodes_marker_enhancers_eqtl_annotation.txt.gz",header = T,sep = "\t") %>% as.data.frame()
nodes_anno <-nodes_anno[,c("Hotspot","eQTL_Source","H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3","CTCF","CHROMATIN_Accessibility","TFBS","Enhancer")]
node_count <-lapply(3:13,function(i){sum(!is.na(nodes_anno[,i]))})
nodes_anno_count =data.frame(markers=colnames(nodes_anno)[3:13],number=unlist(node_count))
nodes_anno_count$fraction= nodes_anno_count$number /nrow(nodes_anno)
write.table(nodes_anno_count,"10_nodes_anno_count.txt",row.names = F, col.names = T,quote =F,sep="\t")

anno <-left_join(hotspot_g,nodes_anno,by="Hotspot")
#-------------------edges
gene_tissue <- read.table("./edges_annotation/07_annotation_hotspot_egene_tissue.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(gene_tissue) <-c("Hotspot","egene","Hotspot_egene_source")
gene_tissue$Hotspot_gene <- paste(gene_tissue$Hotspot,gene_tissue$egene,sep=":")
gene_tissue <-gene_tissue[,-c(1:2)]
enhancer_gene <- read.table("./edges_annotation/success_enhancer_target_anno.bed.gz",header = T,sep = "\t") %>% as.data.frame()
ENH_ENH <- read.table("./edges_annotation/success_ENH_ENH_hotspot_egene.bed.gz",header = T,sep = "\t") %>% as.data.frame()
TSS_ENH <- read.table("./edges_annotation/success_TSS_ENH_hotspot_egene.bed.gz",header = T,sep = "\t") %>% as.data.frame()
TSS_TSS <- read.table("./edges_annotation/success_TSS_TSS_hotspot_egene.bed.gz",header = T,sep = "\t") %>% as.data.frame()
TAD <- read.table("./edges_annotation/09_Adjust_success_egene_hotspot_TAD.bed.gz",header = T,sep = "\t") %>% as.data.frame()
# write.table(anno,"nodes_marker_enhancers_eqtl_annotation.txt",row.names = F, col.names = T,quote =F,sep="\t")
#-------------------interaction
co_ex_interaction <- read.table("02_hotspot_target_gene_reactomeFI_co-expression.bed.gz",header = T,sep = "\t") %>% as.data.frame()

anno$Hotspot_gene <-paste(anno$Hotspot,anno$egene,sep=":")
anno <-left_join(anno,gene_tissue,by="Hotspot_gene")
colnames(enhancer_gene)[1] <-"Hotspot_gene"
anno <-left_join(anno,enhancer_gene,by="Hotspot_gene")
anno <-left_join(anno,ENH_ENH,by="Hotspot_gene")
anno <-left_join(anno,TSS_ENH,by="Hotspot_gene")
anno <-left_join(anno,TSS_TSS,by="Hotspot_gene")
anno <-left_join(anno,TAD,by="Hotspot_gene")

edges_count <-lapply(20:24,function(i){sum(!is.na(anno[,i]))})
edges_anno_count =data.frame(Types=colnames(anno)[20:24],number=unlist(edges_count))
edges_anno_count$fraction= edges_anno_count$number /nrow(anno)
write.table(edges_anno_count,"10_edges_anno_count.txt",row.names = F, col.names = T,quote =F,sep="\t")
length(unique(anno$egene)) #33076

aa <-anno[,c("egene","Enhancer_target","ENH_ENH","TSS_ENH","TSS_TSS","TAD")]
gene<-lapply(2:6,function(i){unique(aa[which(!is.na(aa[,i])),"egene"])})
gene_count<-lapply(2:6,function(i){length(unique(aa[which(!is.na(aa[,i])),"egene"]))})
gene_anno_count = data.frame(Types=colnames(aa)[2:6],number=unlist(gene_count))
write.table(gene_anno_count,"10_genes_anno_count.txt",row.names = F, col.names = T,quote =F,sep="\t")
anno <-anno %>% dplyr::select(-c(Hotspot,Hotspot_gene))
write.table(anno,"10_nodes_edges_annotation.txt",row.names = F, col.names = T,quote =F,sep="\t")