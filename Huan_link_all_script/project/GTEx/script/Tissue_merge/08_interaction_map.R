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

TSS<-read.table("07_TSS_in_win1kb.bed.gz",header = F,sep = "\t") %>% as.data.frame()
eqtl<-read.table("07_eqtl_in_win1kb.bed.gz",header = F,sep = "\t") %>% as.data.frame()

colnames(TSS) <-c("TSS_chr","TSS_start","TSS_end","ENSG","TSS_winnum")
colnames(eqtl) <-c("eQTL_chr","eQTL_start","eQTL_end","ENSG","eQTL_winnum")
inta <-inner_join(TSS,eqtl,by="ENSG")

inta <-unique(inta)


fa <-inta %>%group_by(eQTL_winnum,TSS_winnum)%>%summarise(count =n())%>% as.data.frame()

e <-setdiff(c(1:10000),unique(fa$eQTL_winnum))
S <-setdiff(c(1:10000),unique(fa$TSS_winnum))
ee <-data.frame(eQTL_winnum=setdiff(c(1:10000),unique(fa$eQTL_winnum)),TSS_winnum=1,count=0 )
ss <-data.frame(eQTL_winnum=1,TSS_winnum=setdiff(c(1:10000),unique(fa$TSS_winnum)),count=0 )
ffa <-bind_rows(fa,ee,ss)
x <- dcast(ffa,eQTL_winnum~TSS_winnum)
aaa <-x
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


# color2 = colorRampPalette(c('#FFFFFF','#FE0000'))(50)
color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
pheatmap(r3,cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = F,show_rownames = F,show_colnames = F,file="1MB_win_1kb.png",dpi=120)


#----------

# pdf("071_marker_annotation_ratio_heatmap_filter_3103833.pdf")
# pheatmap(m,cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0,display_numbers = TRUE,show_rownames = F,show_colnames = F) #,cellwidth= 60
# dev.off()


# m[is.na(m)] <-0