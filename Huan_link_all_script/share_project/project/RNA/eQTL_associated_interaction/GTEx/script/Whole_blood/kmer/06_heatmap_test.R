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

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/factor_anno")
All <-read.table("../communities.bed.gz",header = F,sep = "\t") %>% as.data.frame()
CHROMATIN_Accessibility <-read.table("CHROMATIN_Accessibility.bed.gz",header = F,sep = "\t") %>% as.data.frame()
CTCF<-read.table("CTCF.bed.gz",header = F,sep = "\t") %>% as.data.frame()
TFBS<-read.table("TFBS.bed.gz",header = F,sep = "\t") %>% as.data.frame()

H3K27ac<-read.table("H3K27ac.bed.gz",header = F,sep = "\t") %>% as.data.frame() #active
H3K27me3<-read.table("H3K27me3.bed.gz",header = F,sep = "\t") %>% as.data.frame() #repressed
H3K36me3<-read.table("H3K36me3.bed.gz",header = F,sep = "\t") %>% as.data.frame() #transcribed
H3K4me1<-read.table("H3K4me1.bed.gz",header = F,sep = "\t") %>% as.data.frame() #enhancer
H3K4me3<-read.table("H3K4me3.bed.gz",header = F,sep = "\t") %>% as.data.frame() #promoter
H3K9ac<-read.table("H3K9ac.bed.gz",header = F,sep = "\t") %>% as.data.frame() #active
H3K9me3<-read.table("H3K9me3.bed.gz",header = F,sep = "\t") %>% as.data.frame() #gene silencing

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
cluster_number <-c(0,1,2,3,5)
count <-c(20486,7129,4818,4703,830)

# all_anno <-bind_rows(CHROMATIN_Accessibility,CTCF,TFBS,H3K27ac,H3K36me3,H3K4me1,H3K4me3,H3K9ac,H3K27me3,H3K9me3)
all_anno <-bind_rows(CHROMATIN_Accessibility,TFBS,H3K27ac,H3K36me3,H3K4me1,H3K4me3,H3K9ac,H3K27me3,H3K9me3)
colnames(all_anno)[4] <-"cluster_number"
all_anno$pos <-paste(all_anno$V1,all_anno$V2,all_anno$V3,sep="_")
anno_count <-all_anno%>%group_by(cluster_number,Class)%>%summarise(factor=n())%>%as.data.frame()

anno_count_all <-all_anno%>%group_by(Class)%>%summarise(factor=n())%>%as.data.frame()



cluster <-data.frame(cluster_number,count)

final_data <-left_join(anno_count, cluster, by = "cluster_number")
final_data$factor<-as.numeric(as.character(final_data$factor))
final_data$ratio <- final_data$factor/final_data$count
final_data2 <-filter(final_data,cluster_number!=4)

figure_used <-acast(final_data2[,c(1,2,5)],Class~cluster_number)
aaa <-data.frame(figure_used)
aaa<-aaa[c("CHROMATIN Accessibility","TFBS","H3K27ac","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K27me3","H3K9me3"),]
normalization<-function(x){
      return((x -min(x)) / (max(x)-min(x)))}

# normalization<-function(x){
#       return(x /max(x))}

# BBB<-apply(figure_use,1,normalization)%>%data.frame()%>%t()%>%as.matrix()

BBB<-apply(aaa,1,normalization)%>%data.frame()%>%t()
colnames(BBB)<-cluster_number
color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/")
pdf("./figure/06_kmer_heatmap.pdf")
pheatmap(as.matrix(BBB),cluster_rows = FALSE, cluster_cols = FALSE,color= color2,angle_col = 0)
dev.off()



# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/")

# pheatmap(BBB,cluster_rows = FALSE, cluster_cols = FALSE,color= color2)
# dev.off()
# pdf(figure_name)
# title = "qqqq"
# pheatmap(BBB,cluster_rows = FALSE, cluster_cols = FALSE,color = color2,fontsize_row=1.5,fontsize_col=4,cellwidth=7,fontsize=5, main =title)
# dev.off()
