library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(pheatmap)
library(reshape2)

setwd("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/annotation/extend/")
All <-read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()

CHROMATIN_Accessibility <-read.table("CHROMATIN_Accessibility_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()
CTCF<-read.table("CTCF_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()
TFBS<-read.table("TFBS_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()

H3K27ac<-read.table("H3K27ac_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #active
H3K27me3<-read.table("H3K27me3_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #repressed
H3K36me3<-read.table("H3K36me3_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #transcribed
H3K4me1<-read.table("H3K4me1_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #enhancer
H3K4me3<-read.table("H3K4me3_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #promoter
H3K9ac<-read.table("H3K9ac_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #active
H3K9me3<-read.table("H3K9me3_Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #gene silencing

# Active_histone <-bind_rows(H3K27ac,H3K36me3,H3K4me1,H3K4me3,H3K9ac)
# Silent_histone <-bind_rows(H3K27me3,H3K9me3)


CHROMATIN_Accessibility$Class<-"CHROMATIN Accessibility"
CTCF$Class<-"CTCF"
TFBS$Class<-"TFBS"
H3K27ac$Class<-"H3K27ac"
H3K36me3$Class<-"H3K36me3"
H3K4me1$Class<-"H3K4me1"
H3K4me3$Class<-"H3K4me3"
H3K9ac$Class<-"H3K9ac"
H3K27me3$Class<-"H3K27me3"
H3K9me3$Class<-"H3K9me3"

tmp <-bind_rows(CHROMATIN_Accessibility,CTCF,TFBS,H3K27ac,H3K36me3,H3K4me1,H3K4me3,H3K9ac,H3K27me3,H3K9me3)
colnames(tmp)<-c("chr","start","end","chr1","start1","end1","overlap_bp","Factor")
tmp <-tmp%>%select(chr,start,end,Factor)%>%unique()%>%data.frame()
tmp_count <-group_by(tmp,Factor)%>%summarise(count=n())%>%data.frame()
tmp_count$factor_ratio <- tmp_count$count/nrow(All)

# rownames(tmp_count) <-tmp_count$Factor
figure_used <- tmp_count[,3]%>%as.matrix()
rownames(figure_used) <-tmp_count$Factor 
figure_used<-figure_used[c("CHROMATIN Accessibility","TFBS","CTCF","H3K27ac","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K27me3","H3K9me3"),]

setwd("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/figure/")
color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
pdf("071_marker_annotation_ratio_heatmap.pdf")
pheatmap(as.matrix(figure_used),cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = TRUE,show_rownames = T,cellwidth= 60)
dev.off()