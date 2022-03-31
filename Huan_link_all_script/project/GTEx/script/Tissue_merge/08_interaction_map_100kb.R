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
library(parallel)
library(pheatmap)
# library(mclust)
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) #

#---------------------------------------------------------
library(Hmisc)
setwd("/home/huanhuan/project/GTEx/output/Tissue_merge/Cis_eQTL/interaction_heatmap/")


# setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/")
# TSS<-read.table("07_TSS_in_win100.bed.gz",header = F,sep = "\t") %>% as.data.frame()
# eqtl<-read.table("07_eqtl_in_win100.bed.gz",header = F,sep = "\t") %>% as.data.frame()

TSS<-read.table("07_region_half_chr1_hotspot_eqtl_TSS_in_win100000.bed.gz",header = F,sep = "\t") %>% as.data.frame()
eqtl<-read.table("07_region_half_chr1_hotspot_eqtl_in_win100000.bed.gz",header = F,sep = "\t") %>% as.data.frame()
win<-read.table("07_region_half_chr1_win100000.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(TSS) <-c("TSS_chr","TSS_start","TSS_end","ENSG","TSS_winnum")
colnames(eqtl) <-c("eQTL_chr","eQTL_start","eQTL_end","ENSG","eQTL_winnum")
inta <-inner_join(TSS,eqtl,by="ENSG")

inta <-unique(inta)
# inta$win_distance <-inta$TSS_winnum - inta$eQTL_winnum
# inta2 <-filter(inta, abs(inta$win_distance) <=10)

# fa <-inta2 %>%group_by(eQTL_winnum,TSS_winnum)%>%summarise(count =n())%>% as.data.frame()

fa <-inta %>%group_by(eQTL_winnum,TSS_winnum)%>%summarise(count =n())%>% as.data.frame()
L=nrow(win)
ee <-data.frame(eQTL_winnum=setdiff(c(1:L),unique(fa$eQTL_winnum)),TSS_winnum=fa[1,"TSS_winnum"],count=0 )
ss <-data.frame(eQTL_winnum=fa[1,"eQTL_winnum"],TSS_winnum=setdiff(c(1:L),unique(fa$TSS_winnum)),count=0 )
ffa <-bind_rows(fa,ee,ss)
sffa <-ffa[order(-ffa$count),]
# sffa[1,"count"]=950
x <- dcast(sffa,eQTL_winnum~TSS_winnum)
rownames(x) <-x[,1]
x <-x[,-1]
x[is.na(x)] <-0

#-----
r1 <- as.matrix(x)
r1[lower.tri(r1)] <- 0
r1 <- t(r1)

r2 <- as.matrix(x)
r2[upper.tri(r2)] <- 0

r3 = r1 + r2
r3 = r3 + t(r3) 
diag(r3)=diag(r2)
r3[which(r3>0)]=1

# color2 = colorRampPalette(c('#FFFFFF','#FE0000'))(50)
color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
pheatmap(r3,cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="./figure/half_chr1_win_100kb.png",dpi=300)

color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
pheatmap(r3[1:500,1:500],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="./figure/half_chr1_win_100kb_500.png",dpi=300)
pheatmap(r3[1:300,1:300],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="./figure/half_chr1_win_100kb_300.png",dpi=300)
pheatmap(r3[200:500,200:500],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="./figure/half_chr1_win_100kb_200_500.png",dpi=300)
pheatmap(r3[300:600,300:600],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="./figure/half_chr1_win_100kb_300_600.png",dpi=300)
pheatmap(r3[300:400,300:400],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,border=FALSE,file="./figure/half_chr1_win_100kb_300_400.png",dpi=300)
pheatmap(r3[350:400,350:400],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,border=FALSE,file="./figure/half_chr1_win_100kb_350_400.png",dpi=300)
pheatmap(r3[480:520,480:520],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,border=FALSE,file="./figure/half_chr1_win_100kb_480_520.png",dpi=300)
# pheatmap(r3[3500:4000,3500:4000],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="half_chr1_win_10kb_3500_4000.png",dpi=300)
# pheatmap(r3[3500:3750,3500:3750],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="half_chr1_win_10kb_3500_3750.png",dpi=300)

# pheatmap(r3[5000:10000,5000:10000],cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="half_chr1_win_10kb_5000_10000.png",dpi=300)

#----------
# x3 <- 
# a <-ffa[1:40,]
# # a$eQTL_winnum <-paste0("r",a$eQTL_winnum)
# # a$TSS_winnum <-paste0("c",a$TSS_winnum)
# x2 <- dcast(a,eQTL_winnum~TSS_winnum)