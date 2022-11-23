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

setwd("/home/huanhuan/project/eQTL_Catalogue/output/annotation/extend_18snp/")
All <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()

CHROMATIN_Accessibility <-read.table("CHROMATIN_Accessibility_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()
CTCF<-read.table("CTCF_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()
TFBS<-read.table("TFBS_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()

H3K27ac<-read.table("H3K27ac_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #active
H3K27me3<-read.table("H3K27me3_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #repressed
H3K36me3<-read.table("H3K36me3_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #transcribed
H3K4me1<-read.table("H3K4me1_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #enhancer
H3K4me3<-read.table("H3K4me3_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #promoter
H3K9ac<-read.table("H3K9ac_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #active
H3K9me3<-read.table("H3K9me3_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame() #gene silencing

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
colnames(tmp)<-c("chr","start","end","chr1","start1","end1","overlap","Factor")
tmp$hotspot <- paste(tmp$chr,tmp$start,tmp$end,sep="_")
tmp$overlap <-1
tmp <-tmp%>%select(hotspot,Factor,overlap)%>%unique()%>%data.frame()
#==================================================================================
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/markers")
all_hotspot <- read.table("../GWAS/17_whole_genome_leiden_pca3_k50_resolution1e-04_sorted.bed.gz",header = F,sep = "\t") %>% as.data.frame()
all_hotspot$hotspot <-paste(all_hotspot$V1,all_hotspot$V2,all_hotspot$V3,sep="_")
colnames(all_hotspot)[4] <-"cluster"
all_hotspot <-all_hotspot[,c("hotspot","cluster")]
dat <- left_join(all_hotspot,tmp,by="hotspot")
fdat<-dcast(dat,hotspot+cluster~Factor)
fdat <-fdat[,-13]
fdat[is.na(fdat)]<-0

library(UpSetR)
pdf("28_upset_plot_marker_all_hotspot.pdf",width=13,height=7)
upset(fdat,mainbar.y.label="Number of hotspot",sets.x.label="Number of hotspot",sets = c("CHROMATIN Accessibility","TFBS","CTCF","H3K27ac","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K27me3","H3K9me3"))
dev.off()


ffdat <-head(fdat,1000)
library(UpSetR)
pdf("28_upset_plot_marker_all_hotspot1.pdf",width=13,height=7)
upset(ffdat,mainbar.y.label="Number of hotspot",sets.x.label="Number of hotspot",sets = c("H3K4me3","H3K9ac","H3K27me3","H3K9me3"))
dev.off()


aa <-filter(fdat,H3K4me3==1&H3K9ac==0&H3K27me3==0&H3K27ac==0&H3K9me3==0&H3K36me3==0&H3K4me1==0)

tmp_count <-group_by(tmp,Factor)%>%summarise(count=n())%>%data.frame()
tmp_count$factor_ratio <- tmp_count$count/nrow(All)

# rownames(tmp_count) <-tmp_count$Factor
figure_used <- tmp_count[,3]%>%as.matrix()
rownames(figure_used) <-tmp_count$Factor 
figure_used<-figure_used[c("CHROMATIN Accessibility","TFBS","CTCF","H3K27ac","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K27me3","H3K9me3"),]

setwd("/home/huanhuan/project/eQTL_Catalogue/script/figure/")
color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
pdf("071_marker_annotation_ratio_heatmap.pdf")
pheatmap(as.matrix(figure_used),cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = TRUE,show_rownames = T,cellwidth= 60)
dev.off()