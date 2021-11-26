
source("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/script/OmicCircos_circos_refine_revise_color.R")
library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(tidyverse)
library(reshape2)
library(Seurat)
library("OmicCircos")
data("UCSC.hg19.chr")


# cancers <-c("BRCA","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","STAD","THCA")
# cancers <-c("BRCA","KIRC","KIRP","LAML","LIHC","LUAD","OV","PAAD","PRAD","STAD","THCA")
tissues <-  c("Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Prostate","Thyroid","Whole_Blood")
re_file<-function(tissue=NULL){
    path = paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/hotspot_cis_eQTL/interval_18")
    setwd(path)
    tissue2 <-str_replace(tissue,"Whole_Blood","whole_blood")
    # read_file = paste0(tissue2,"_segment_hotspot_cutoff_0.176_makewin_1MB.bed.gz")
    # read_file = paste0(tissue2,"_segment_hotspot_cutoff_0.176_makewin_20KB.bed.gz")
    # read_file = paste0(tissue2,"_segment_hotspot_cutoff_0.176_makewin_10MB.bed.gz")
    # read_file = paste0(tissue2,"_segment_hotspot_cutoff_0.176_makewin_1KB.bed.gz")
    # read_file = paste0(tissue2,"_segment_hotspot_cutoff_0.176_makewin_10KB.bed.gz")
    # read_file = paste0(tissue2,"_segment_hotspot_cutoff_0.176_makewin_5KB.bed.gz")
    read_file = paste0(tissue2,"_segment_hotspot_cutoff_0.176_makewin_261BP.bed.gz")
    org <- read.table(read_file,header = F,sep = "\t") %>% as.data.frame()
    colnames(org) <-c("chr","start","end","overlap")
    org$tissue <-tissue
    return(org)
}

tmp <-lapply(tissues,re_file)
rs <-do.call(rbind,tmp)

sum_count <-rs%>%group_by(chr,start,end,tissue)%>%summarize(hot_sum=sum(overlap))%>%data.frame()
sum_count$ratio <-sum_count$hot_sum/(sum_count$end - sum_count$start)
sum_count$pos <- paste(sum_count$chr,sum_count$start,sum_count$end,sep="_")
# sum_count <- sum_count%>%dplyr::select(-hot_sum)
sum_count <- sum_count[,c(4,6,7)]
sum_count <-sum_count[,c("pos","tissue","ratio")]
final_sum_count <- acast(sum_count,pos~tissue)%>%data.frame()

final_sum_count$pos <-rownames(final_sum_count)

tmp2 <- str_split_fixed(rownames(final_sum_count ),pattern='_',n=3)

final_sum_count <-add_column(final_sum_count, chr = as.data.frame(tmp2)[,1], .before = 1)
final_sum_count <-add_column(final_sum_count, start = as.data.frame(tmp2)[,2], .before = 2)
final_sum_count <-add_column(final_sum_count, end = as.data.frame(tmp2)[,3], .before = 3)

final_sum_count <-final_sum_count%>%dplyr::select(-pos)
rownames(final_sum_count) <-NULL
final_sum_count[is.na(final_sum_count)]=0
#------------------
final_sum_count1 <-as.matrix(final_sum_count)



library(circlize)
circos.clear()
pdf("16_circos_hotspot_densityaaaa.pdf")
par(mar = c(2, 2, 2, 2))
circos.par(start.degree = 90)
circos.initializeWithIdeogram(species= "hg19",chromosome.index = paste0("chr", 1:22))
aaa <-head(final_sum_count,10)
circos.lines(x=aaa$start,y=aaa$Breast_Mammary_Tissue,sector.index="chr10")




final_sum_count2 <-final_sum_count%>%dplyr::select(-end) %>%as.matrix()
dat <- UCSC.hg19.chr
#dat$chrom <- gsub("chr", "",dat$chrom)
indx = which(dat$chrom != "chrX" & dat$chrom != "chrY")
seg.d = dat[indx, ]
seg.name = paste("chr", seq(1, 22, by = 1), sep = "")
seg.c = segAnglePo(seg.d, seg = seg.name)
# colors <- rainbow(10, alpha=0.5);
# colors <- rainbow(22, alpha = 0.8)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_total/")
# pdf(file = "circos_heatmap_split_1MB.pdf", height = 8, width = 8, compress = TRUE)
# pdf(file = "circos_heatmap_split_10MB.pdf", height = 8, width = 8, compress = TRUE)
# pdf(file = "circos_heatmap_split_20KB.pdf", height = 8, width = 8, compress = TRUE)
# pdf(file = "circos_heatmap_split_1KB.pdf", height = 8, width = 8, compress = TRUE)
# pdf(file = "circos_heatmap_split_10KB.pdf", height = 8, width = 8, compress = TRUE)
# pdf(file = "circos_heatmap_split_5KB.pdf", height = 8, width = 8, compress = TRUE)
pdf(file = "circos_heatmap_split_261BP.pdf", height = 8, width = 8, compress = TRUE)
par(mar= c(1,1,1,1)) #设置画图参数
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "")
# circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos_refine(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE,col ='#2171b5')
color2<-c('#c6dbef','#6baed6','#2171b5')
# circos(R = 20, cir = seg.c, type= "heatmap", W = 380, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE,col ="red")
circos_refine(R = 100, cir = seg.c, type= "heatmap", W = 280, mapping = final_sum_count1,B = T, huan_col=color2,col.v = 4,col.bar=T)
dev.off()

#-----------------
Breast_Mammary_Tissue <- final_sum_count1[,c(1:4)]
pdf(file = "circos_heatmap_split_261BP_Breast_Mammary_Tissue.pdf", height = 8, width = 8, compress = TRUE)
par(mar= c(1,1,1,1)) #设置画图参数
plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "")
# circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
circos_refine(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE,col ='#2171b5')
color2<-c('#c6dbef','#6baed6','#2171b5')
# circos(R = 20, cir = seg.c, type= "heatmap", W = 380, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE,col ="red")
circos_refine(R = 100, cir = seg.c, type= "heatmap", W = 280, mapping = Breast_Mammary_Tissue ,B = T, huan_col=color2,col.v = 4,col.bar=T)
dev.off()

# par(mar= c(1,1,1,1)) #设置画图参数
# p1 <-plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "")
# # circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
# p1<- p1+circos_refine(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE,col ='#2171b5')
# color2<-c('#c6dbef','#6baed6','#2171b5')
# # circos(R = 20, cir = seg.c, type= "heatmap", W = 380, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE,col ="red")
# p1 <- p1+circos_refine(R = 100, cir = seg.c, type= "heatmap", W = 280, mapping = final_sum_count1,B = T, huan_col=color2,col.v = 4,col.bar=T)



# ggsave("circos_heatmap_split_261BP.png", p1,dpi=300)
# ggsave(paste0(savePATH,"ini_mut_count.png"),p1,dpi=300)
#------------------
# final_sum_count3 <-final_sum_count[,c(1:4)] %>%as.matrix()

# pdf(file = "aaa.pdf", height = 8, width = 8, compress = TRUE)
# # pdf(file = "circos_heatmap_split_20KB.pdf", height = 8, width = 8, compress = TRUE)
# par(mar= c(1,1,1,1)) #设置画图参数
# plot(c(1,800), c(1,800), type= "n", axes = FALSE, xlab = "", ylab = "")
# # circos(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, col = colors, scale = TRUE, print.chr.lab = TRUE, cex = 2)
# circos_refine(R = 400, cir = seg.c, mapping = seg.d, type = "chr", W = 5, scale = TRUE, print.chr.lab = TRUE,col ='#2171b5')
# color2<-c('#c6dbef','#6baed6','#2171b5')
# # circos(R = 20, cir = seg.c, type= "heatmap", W = 380, mapping = data1, lwd = 0.01, col.v = 4, cluster=FALSE, col.bar=TRUE,col ="red")
# circos_refine(R = 100, cir = seg.c, type= "heatmap", W = 280, mapping = final_sum_count3,B = T, huan_col=color2,col.v = 4,col.bar=T)
# dev.off()

