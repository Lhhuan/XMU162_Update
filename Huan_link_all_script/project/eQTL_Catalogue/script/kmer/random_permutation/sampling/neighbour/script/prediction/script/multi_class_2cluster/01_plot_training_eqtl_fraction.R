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
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/multi_class_2cluster/output/")
posi <- fread("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_0.05_pvalue.bed.gz",header = F,sep = "\t") %>% as.data.frame()
all_posi <-  read.csv("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/GC_content/gc_content.bed.gz",header = T,sep = "\t") %>% as.data.frame()

colnames(all_posi) <- c("hchr","hstart","h_end","gc_content","class")

colnames(posi) <- c("e_chr","e_start","e_end","egene","pvalue","hchr","hstart","h_end")
posi1 <-posi%>%group_by(e_chr,e_start,e_end,hchr,hstart,h_end)%>%summarise(mim_p=min(pvalue))%>%data.frame()
posi2 <-left_join(posi1,all_posi[,c(1:3,5)],by=c("hchr","hstart","h_end"))

nega <- fread("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/01_all_negative_egene_0.05_pvalue.bed.gz",header = F,sep = "\t") %>% as.data.frame()
all_nega <-fread("/share/data0/QTLbase/huan/eQTL_Catalogue/original_random/extend_18snp/emp0.176/0/1_resemble_Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t")%>% as.data.frame()
colnames(nega) <- c("e_chr","e_start","e_end","egene","pvalue","hchr","hstart","h_end")
nega1 <-nega%>%group_by(e_chr,e_start,e_end,hchr,hstart,h_end)%>%summarise(mim_p=min(pvalue))%>%data.frame()
nega1$class <- 0
all_eqtl <-bind_rows(posi2,nega1)
all_nega$class <-0
colnames(all_nega) <- c("hchr","hstart","h_end","class")

all_dat <- bind_rows(all_nega,all_posi[,c("hchr","hstart","h_end","class")])
all_dat_n <- data.frame(table(all_dat$class))


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
    dat1 <-filter(all_eqtl,mim_p <cutoff)
    dat2 <- unique(dat1[,c(4:6,8)])
    dat3 <- data.frame(table(dat2$class))
    org5 <- left_join(dat3,all_dat_n,by="Var1")
    colnames(org5) <-c("class","N_cutoff","N_total")
    org5$class <-as.character(org5$class)
    org5[4,] <- data.frame("1_2",org5[2,2]+org5[3,2],org5[2,3]+org5[3,3])
    org5$cutoff <- cutoff
    rs <-bind_rows(rs,org5)
    print(cutoff)
}
rs <-rs[-23,]

rs$class <- paste0("C",rs$class)
rs$cutoff <-as.factor(rs$cutoff)
rs$prop <- rs$N_cutoff/rs$N_total
rs1 <-rs
rs[23,5] <-1
rs$class <-factor(rs$class,levels=c("C0","C1_2","C1","C2"))

p1 <- ggplot(rs,mapping=aes(x=cutoff,y=prop,fill=class))+geom_bar(stat = 'identity',position=position_dodge(0.9))+ #

# scale_fill_manual(values=c("#A593E0","#F68657","#F6B352","#84B1ED"))+
scale_fill_manual(values=c("#A593E0","#F68657","#84B1ED","#F6B352"))+
p_theme+
labs(x="",y="proportion of segments",fill="")
pdf("01_eqtl_segment_proportion.pdf",height=4,width=5)
print(p1)
dev.off()