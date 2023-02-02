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
# library(mclust)
library(umap)
library(uwot)
library(kernlab)
library(data.table)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
posi <- fread("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05_pvalue.bed.gz",header = F,sep = "\t") %>% as.data.frame()
all_posi <-  fread("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()


colnames(posi) <- c("e_chr","e_start","e_end","egene","pvalue","hchr","hstart","h_end")
posi1 <-posi%>%group_by(e_chr,e_start,e_end,hchr,hstart,h_end)%>%summarise(mim_p=min(pvalue))%>%data.frame()
# posi <-posi[,c("e_chr","e_start","e_end","hchr","hstart","h_end")]%>%unique()

nega <- fread("01_all_negative_egene_0.05_pvalue.bed.gz",header = F,sep = "\t") %>% as.data.frame()
all_nega <-fread("/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t")%>% as.data.frame()
colnames(nega) <- c("e_chr","e_start","e_end","egene","pvalue","hchr","hstart","h_end")
nega1 <-nega%>%group_by(e_chr,e_start,e_end,hchr,hstart,h_end)%>%summarise(mim_p=min(pvalue))%>%data.frame()
#==========================

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))


#===========
# cutoff=1e-3
rs <-data.frame()
# filter_cutoff <- function(){
for (cutoff in c(0.05,1e-2, 1e-3,1e-4,1e-5,5e-8)){
    posi2 <- filter(posi1,mim_p<cutoff)
    posi_n <- unique(posi2[,4:6])
    posi2$class <- "True"
    #=====================
    nega2 <- filter(nega1,mim_p<cutoff)
    nega_n <-unique(nega2[,4:6])
    cover_eqtl=c(nrow(posi_n),nrow(nega_n))
    total=c(nrow(all_posi),nrow(all_nega))
    class=c("True","False")
    prop_t <- data.frame(cover_eqtl=c(nrow(posi_n),nrow(nega_n)),total=c(nrow(all_posi),nrow(all_nega)),class=c("True","False"))
    prop_t$prop <- prop_t$cover_eqtl/prop_t$total
    prop_t$cutoff <- cutoff 
    rs <-bind_rows(rs,prop_t)
    print(cutoff)
}

rs$class <-factor(rs$class,levels=c("True","False"))
rs$cutoff <-as.factor(rs$cutoff)
p1 <- ggplot(rs,mapping=aes(x=cutoff,y=prop,fill=class))+geom_bar(stat = 'identity',position=position_dodge(0.9))+ #
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    scale_fill_manual(values=c("#A593E0","#F68657"))+
    p_theme+
    labs(x="",y="proportion of segments",fill="")
pdf("./figure/013_training_segment_proportion.pdf",height=4,width=5)
print(p1)
dev.off()